test_that("local polynomial estimators retain their reference values", {
  X <- c(seq(-1, -0.1, by = 0.1), seq(0.1, 1, by = 0.1))
  Y <- 1 + 1.5 * (X >= 0) + 0.7 * X - 0.4 * X^2 +
    c(
      -0.20, 0.10, 0.05, -0.10, 0.15, -0.05, 0.08, -0.12, 0.03, 0.07,
      -0.04, 0.11, -0.09, 0.06, -0.02, 0.13, -0.07, 0.04, -0.11, 0.09
    )
  GG <- diag(length(X))

  conventional <- rddinterf:::estimate_conv(0.8, X, Y, GG)
  robust <- rddinterf:::estimate_bc(0.8, 1, X, Y, GG)

  expect_equal(conventional$est, 1.45357142857143, tolerance = 1e-12)
  expect_equal(conventional$var, 0.00582920814244062, tolerance = 1e-12)
  expect_equal(c(conventional$N_l, conventional$N_r), c(7, 7))

  expect_equal(robust$est, 1.37467532467533, tolerance = 1e-12)
  expect_equal(as.numeric(robust$var), 0.00631396069324256, tolerance = 1e-12)
  expect_equal(c(robust$N_l, robust$N_r), c(7, 7))
  expect_equal(
    as.numeric(robust$beta_plus_bc - robust$beta_minus_bc),
    robust$est,
    tolerance = 1e-12
  )
})

test_that("estimators handle subsets and missing values consistently", {
  X <- seq(-1, 1, length.out = 30)
  Y <- X + (X >= 0) + sin(seq_along(X)) / 10
  GG <- diag(length(X))
  keep <- rep(c(TRUE, TRUE, FALSE), 10)
  Y[1] <- NA_real_

  fit <- rddinterf:::estimate_conv(0.9, X, Y, GG, subset = keep)

  expect_true(is.finite(fit$est))
  expect_true(is.finite(fit$var))
  expect_gt(fit$N_l, 0)
  expect_gt(fit$N_r, 0)
})

test_that("the public interface returns all requested effects", {
  set.seed(123)
  n <- 80
  Xi <- runif(n, -1, 1)
  neighbors <- lapply(seq_len(n), function(i) {
    c(if (i == 1) n else i - 1L, if (i == n) 1L else i + 1L)
  })
  W <- diag(n)
  for (i in seq_len(n)) W[i, neighbors[[i]]] <- 1
  W <- (W + t(W) > 0) * 1
  Gi <- rddinterf:::one_treated(as.integer(Xi >= 0), neighbors)
  Yi <- 1 + (Xi >= 0) + 0.3 * Gi + 0.5 * Xi + rnorm(n, sd = 0.2)
  contrasts <- list(
    direct = list(c(1, 0), c(0, 0)),
    spillover = list(c(0, 1), c(0, 0))
  )

  fit <- rddinterf(
    Yi, Xi, neighbors, W, contrasts, "one.treated",
    estimator = "conventional", h = 0.8, b = 1
  )

  expect_s3_class(fit, "rddinterf")
  expect_named(
    fit$estimates,
    c("ov_direct", "ov_indirect", "direct", "spillover")
  )
  expect_identical(fit$estimator, "conventional")
  expect_equal(dim(summary(fit)), c(4, 7))
  expect_true(all(vapply(fit$estimates, function(x) {
    isTRUE(all.equal(x$mean_right - x$mean_left, x$est, tolerance = 1e-12))
  }, logical(1))))

  printed <- capture.output(returned <- print(fit))
  expect_identical(returned, fit)
  expect_true(any(grepl("Conventional local polynomial", printed, fixed = TRUE)))
  expect_true(any(grepl("Exposure mapping: one.treated", printed, fixed = TRUE)))
  expect_true(any(grepl("spillover", printed, fixed = TRUE)))
  expect_false(any(grepl("dstar|boundary|Yi", printed)))
})

test_that("the public interface reports malformed inputs clearly", {
  Yi <- 1:4
  Xi <- c(-1, -0.5, 0.5, 1)
  neighbors <- list(2L, 1L, 4L, 3L)
  W <- diag(4)
  contrast <- list(effect = list(c(1, 0), c(0, 0)))

  expect_error(
    rddinterf(Yi, Xi, neighbors, W, unname(contrast), "one.treated"),
    "named list"
  )
  expect_error(
    rddinterf(Yi, Xi, neighbors[-1], W, contrast, "one.treated"),
    "same length"
  )
  expect_error(
    rddinterf(Yi, Xi, neighbors, diag(3), contrast, "one.treated"),
    "square matrix"
  )
  expect_error(
    rddinterf(Yi, c(NA, Xi[-1]), neighbors, W, contrast, "one.treated"),
    "Xi must not"
  )
})

test_that("pooled contrasts warn when boundary codimension varies", {
  Xi <- c(-0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8)
  Yi <- seq_along(Xi)
  neighbors <- list(5L, 1L, c(5L, 6L), c(1L, 2L), 1L, 2L, c(1L, 2L), c(5L, 6L))
  W <- diag(length(Xi))
  contrast <- list(indirect = list(c(0, 1), c(0, 0)))

  expect_warning(
    suppressMessages(rddinterf(
      Yi, Xi, neighbors, W, contrast, "frac.treated",
      h = 1, b = 1,
      overall_direct = FALSE, overall_indirect = FALSE
    )),
    "minimum-codimension strata"
  )

  expect_no_warning(
    suppressMessages(rddinterf(
      Yi, Xi, neighbors, W, contrast, "one.treated",
      h = 1, b = 1,
      overall_direct = FALSE, overall_indirect = FALSE
    ))
  )
})
