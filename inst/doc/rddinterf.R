## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(rddinterf)

set.seed(123)
n <- 80
Xi <- runif(n, -1, 1)
neighbors <- lapply(seq_len(n), function(i) {
  c(if (i == 1) n else i - 1L, if (i == n) 1L else i + 1L)
})

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

