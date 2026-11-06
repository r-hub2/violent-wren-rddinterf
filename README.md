# rddinterf

`rddinterf` estimates causal effects in sharp regression discontinuity designs
when units may interfere through a network. It implements the distance-based
local polynomial estimators and network-robust inference developed by Elena
Dal Torrione, Tiziano Arduini, and Laura Forastiere in *Regression
Discontinuity Designs Under Interference*.

The package has two main functions:

- `get_boundary()` constructs an effective-treatment boundary.
- `rddinterf()` estimates overall and contrast-specific effects.

## Example

```r
library(rddinterf)

set.seed(123)
n <- 80
Xi <- runif(n, -1, 1)

# A ring network: every unit has two neighbors.
neighbors <- lapply(seq_len(n), function(i) {
  c(if (i == 1) n else i - 1L, if (i == n) 1L else i + 1L)
})

# Dependency graph, including self-dependence on the diagonal.
W <- diag(n)
for (i in seq_len(n)) W[i, neighbors[[i]]] <- 1
W <- (W + t(W) > 0) * 1

Di <- as.integer(Xi >= 0)
Gi <- as.numeric(vapply(neighbors, function(j) any(Di[j] == 1), logical(1)))
Yi <- 1 + Di + 0.3 * Gi + 0.5 * Xi + rnorm(n, sd = 0.2)

contrasts <- list(
  direct = list(c(1, 0), c(0, 0)),
  spillover = list(c(0, 1), c(0, 0))
)

fit <- rddinterf(
  Yi, Xi, neighbors, W, contrasts,
  exposure = "one.treated",
  estimator = "conventional",
  h = 0.8,
  b = 1
)

fit
```

Printing the fitted object gives a compact results table. Use `summary(fit)` to
obtain the underlying numeric summary matrix.

The supported exposure mappings are `"one.treated"`, `"sum.treated"`,
`"frac.treated"`, and `"majority"`. See `vignette("rddinterf")` for the input
conventions and a fuller explanation of the returned estimates.

When a contrast pools neighborhood strata with different boundary
codimensions, `rddinterf()` warns that the pooled estimator asymptotically targets
the minimum-codimension strata. Estimate each neighborhood stratum separately
when effects for every stratum are desired.

## Reference

Dal Torrione, E., Arduini, T., and Forastiere, L. (2024). “Regression
Discontinuity Designs Under Interference.” arXiv:2410.02727.
