#' Estimate regression discontinuity effects under interference
#'
#' Estimates overall direct effects, overall indirect effects, and
#' specific interference effects in regression discontinuity designs with
#' interference using distance-based local polynomial estimators and
#' network-robust variance estimators.
#'
#' The function supports two estimators:
#'
#'   \describe{
#'     \item{`"conventional"`}{Local polynomial estimator without bias correction.}
#'     \item{`"robust"`}{Bias-corrected local polynomial estimator with robust inference.}
#'   }
#'
#' The own treatment indicator is computed internally from `Xi` and `cutoff`.
#' The neighbors' running-variable lists and neighborhood treatments are also
#' computed internally from `Xi`, `neighbors`, and `exposure`.
#'
#' If a contrast pools active neighborhood strata with different boundary
#' codimensions, the function issues a warning. In that setting, the pooled
#' estimator converges to an effect over the minimum-codimension strata. Effects
#' for every stratum should instead be estimated separately and aggregated with
#' the appropriate density weights.
#'
#' @examples
#' set.seed(123)
#' n <- 80
#' Xi <- runif(n, -1, 1)
#' neighbors <- lapply(seq_len(n), function(i) {
#'   c(if (i == 1) n else i - 1L, if (i == n) 1L else i + 1L)
#' })
#' W <- diag(n)
#' for (i in seq_len(n)) W[i, neighbors[[i]]] <- 1
#' W <- (W + t(W) > 0) * 1
#' Di <- as.integer(Xi >= 0)
#' Gi <- as.numeric(vapply(neighbors, function(j) any(Di[j] == 1), logical(1)))
#' Yi <- 1 + Di + 0.3 * Gi + 0.5 * Xi + rnorm(n, sd = 0.2)
#' contrasts <- list(
#'   direct = list(c(1, 0), c(0, 0)),
#'   spillover = list(c(0, 1), c(0, 0))
#' )
#'
#' fit <- rddinterf(
#'   Yi, Xi, neighbors, W, contrasts, exposure = "one.treated",
#'   estimator = "conventional", h = 0.8, b = 1
#' )
#' fit
#'
#' @param Yi Numeric vector of observed outcomes.
#' @param Xi Numeric vector of individual running variables.
#' @param neighbors List whose `i`th element contains the observation indices
#'   of unit `i`'s neighbors.
#' @param W Binary dependency graph, with rows and columns ordered like `Yi`.
#'   Base and sparse matrices are supported.
#' @param contrasts Named list of effective treatment contrasts to estimate.
#'   Each element should define a treated and control effective treatment.
#' @param exposure Character string specifying the exposure mapping. Supported
#'   values currently include `"one.treated"`, `"sum.treated"`,
#'   `"frac.treated"`, and `"majority"`.
#' @param boundary Optional named list of boundary objects returned by
#'   [get_boundary()]. If `NULL`, boundaries are constructed internally.
#' @param cutoff Numeric cutoff defining treatment assignment. Default is `0`.
#' @param subset Optional logical vector indicating the subset of observations
#'   to use.
#' @param estimator Character string specifying the estimator. Must be one of
#'   `"conventional"` or `"robust"`.
#' @param p Integer specifying the order of the local polynomial estimator.
#' @param q Integer specifying the order of the local polynomial estimator used
#'   for bias correction.
#' @param h Optional manual bandwidth for estimation.
#'   Can be a scalar, named vector, or named list.
#' @param b Optional manual bandwidth for bias correction.
#'   Can be a scalar, named vector, or named list.
#' @param overall_direct Logical; if `TRUE`, compute the overall direct effect.
#' @param overall_indirect Logical; if `TRUE`, compute the overall indirect effect.
#'
#' @returns An object of class `"rddinterf"` with components:
#'
#'   \describe{
#'     \item{`estimates`}{Named list of estimation results.}
#'     \item{`boundary`}{List of boundary objects used for estimation.}
#'     \item{`exposure`}{Exposure mapping used.}
#'     \item{`estimator`}{Estimation method used.}
#'     \item{`cutoff`}{Cutoff used to define treatment.}
#'     \item{`p`}{Order of the main local polynomial fit.}
#'     \item{`q`}{Order of the bias-correction fit.}
#'     \item{`call`}{Matched function call.}
#'   }
#'
#' Each element of `estimates` is either `NA` or a list containing:
#'
#'   \describe{
#'     \item{`est`}{Point estimate.}
#'     \item{`se`}{Estimated standard error.}
#'     \item{`ci`}{95% confidence interval.}
#'     \item{`mean_left`}{Estimated conditional mean at the boundary from the left.}
#'     \item{`mean_right`}{Estimated conditional mean at the boundary from the right.}
#'     \item{`boundary`}{Boundary object or boundary description.}
#'     \item{`dstar`}{Signed running variable used for estimation.}
#'     \item{`Yi`}{Outcome vector in the estimation sample.}
#'     \item{`N_r`}{Effective sample size to the right of the cutoff.}
#'     \item{`N_l`}{Effective sample size to the left of the cutoff.}
#'     \item{`h`}{Bandwidth used for estimation.}
#'     \item{`b`}{Bandwidth used for bias correction.}
#'     \item{`tr`}{Logical indicator for the estimation sample.}
#'     \item{`p`}{Order of the main local polynomial fit.}
#'     \item{`q`}{Order of the bias-correction fit.}
#'   }
#'
#' @export
rddinterf <- function(
    Yi, Xi, neighbors, W,
    contrasts, exposure,
    boundary = NULL,
    cutoff = 0,
    subset = NULL,
    estimator = c("conventional", "robust"),
    p = 1,
    q = p + 1,
    h = NULL,
    b = NULL,
    overall_direct = TRUE,
    overall_indirect = TRUE
) {
  estimator <- match.arg(estimator)

  if (!is.numeric(Yi) || !is.numeric(Xi)) {
    stop("Yi and Xi must be numeric vectors.")
  }

  if (anyNA(Xi) || any(!is.finite(Xi))) {
    stop("Xi must not contain missing or non-finite values.")
  }

  if (!is.list(neighbors)) {
    stop("neighbors must be a list of integer index vectors.")
  }

  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff)) {
    stop("cutoff must be one finite number.")
  }

  if (!is.numeric(p) || length(p) != 1L || p < 0 || p != as.integer(p)) {
    stop("p must be a nonnegative integer.")
  }

  if (!is.numeric(q) || length(q) != 1L || q != as.integer(q)) {
    stop("q must be an integer.")
  }

  if (q < p + 1) {
    stop("q must be at least p + 1.")
  }

  n <- length(Yi)

  if (length(Xi) != n) {
    stop("Xi and Yi must have the same length.")
  }

  if (length(neighbors) != n) {
    stop("neighbors must have the same length as Yi.")
  }

  valid_neighbors <- vapply(neighbors, function(x) {
    is.numeric(x) && all(is.finite(x)) && all(x == as.integer(x)) &&
      all(x >= 1L & x <= n)
  }, logical(1))
  if (!all(valid_neighbors)) {
    stop("Each element of neighbors must contain valid observation indices.")
  }

  if (is.null(dim(W)) || length(dim(W)) != 2L ||
      any(dim(W) != c(n, n))) {
    stop("W must be a square matrix with one row and column per observation.")
  }

  if (!is.null(subset) &&
      (!is.logical(subset) || length(subset) != n || anyNA(subset))) {
    stop("subset must be a logical vector without missing values and have the same length as Yi.")
  }

  if (!is.list(contrasts) || length(contrasts) == 0L ||
      is.null(names(contrasts)) || any(names(contrasts) == "")) {
    stop("contrasts must be a nonempty named list of treatment contrasts.")
  }

  valid_contrasts <- vapply(contrasts, function(x) {
    is.list(x) && length(x) == 2L &&
      all(vapply(x, function(z) is.numeric(z) && length(z) == 2L,
                 logical(1)))
  }, logical(1))
  if (!all(valid_contrasts)) {
    stop("Each treatment contrast must be list(c(d, g), c(dprime, gprime)).")
  }

  if (!is.null(boundary) && length(boundary) != length(contrasts)) {
    stop("boundary and contrasts must have the same length.")
  }

  Di <- as.integer(Xi >= cutoff)
  Xi_centered <- Xi - cutoff
  Xj <- lapply(neighbors, function(idx) Xi_centered[idx])

  Gi <- switch(
    exposure,
    "one.treated" = one_treated(Di, neighbors),
    "sum.treated" = sum_treated(Di, neighbors),
    "frac.treated" = frac_treated(Di, neighbors),
    "majority" = majority(Di, neighbors),
    stop("Unsupported exposure.")
  )

  if (is.null(boundary)) {
    boundary <- lapply(contrasts, function(tr) {
      get_boundary(neighbors, tr, exposure)
    })
    names(boundary) <- names(contrasts)
  }

  Xmin <- vapply(Xj, function(x) {
    if (length(x) == 0L) return(NA_real_)
    x[which.min(abs(x))]
  }, numeric(1))

  estimate <- list()

  if (overall_direct) {
    keep <- if (is.null(subset)) rep(TRUE, n) else subset

    estimate$ov_direct <- fit_one(
      name = "ov_direct",
      Xrun = Xi_centered,
      keep = keep,
      boundary_obj = "Xi = cutoff",
      sbar = 1,
      Y = Yi,
      GG = W,
      estimator = estimator,
      p = p,
      q = q,
      h_manual = h,
      b_manual = b
    )
  }

  if (overall_indirect) {
    keep <- if (is.null(subset)) !is.na(Xmin) else (!is.na(Xmin) & subset)

    estimate$ov_indirect <- fit_one(
      name = "ov_indirect",
      Xrun = Xmin,
      keep = keep,
      boundary_obj = "Xmin = 0",
      sbar = 1,
      Y = Yi,
      GG = W,
      estimator = estimator,
      p = p,
      q = q,
      h_manual = h,
      b_manual = b
    )
  }

  for (k in seq_along(contrasts)) {
    d <- contrasts[[k]][[1]][1]
    g <- contrasts[[k]][[1]][2]
    dprime <- contrasts[[k]][[2]][1]
    gprime <- contrasts[[k]][[2]][2]
    nm <- names(contrasts)[k]

    tr <- (Di == d & Gi == g)
    cr <- (Di == dprime & Gi == gprime)
    keep <- if (is.null(subset)) (tr | cr) else ((tr | cr) & subset)

    bb <- boundary[[k]]
    boundary_sizes <- vapply(bb, function(x) x[[1]]$dim, numeric(1))
    active_boundaries <- boundary_sizes %in%
      unique(vapply(neighbors[keep], length, integer(1)))
    codimensions <- vapply(bb, function(x) x[[1]]$ss, numeric(1))
    if (any(active_boundaries)) codimensions <- codimensions[active_boundaries]
    finite_codimensions <- sort(unique(codimensions[is.finite(codimensions)]))

    if (length(finite_codimensions) > 1L) {
      warning(
        "Contrast '", nm, "' pools neighborhood strata with boundary ",
        "codimensions ", paste(finite_codimensions, collapse = ", "), ". ",
        "When codimension varies, the pooled estimator converges to an effect ",
        "over the minimum-codimension strata.",
        call. = FALSE
      )
    }

    sbar <- min(codimensions)

    dstar <- compute_dstar(
      Xi = Xi_centered,
      Xj = Xj,
      keep = keep,
      d = d,
      g = g,
      dprime = dprime,
      gprime = gprime,
      exposure = exposure
    )

    dstar[cr & keep] <- -abs(dstar[cr & keep])
    dstar[tr & keep] <- abs(dstar[tr & keep])

    estimate[[nm]] <- fit_one(
      name = nm,
      Xrun = dstar,
      keep = keep,
      boundary_obj = bb,
      sbar = sbar,
      Y = Yi,
      GG = W,
      estimator = estimator,
      p = p,
      q = q,
      h_manual = h,
      b_manual = b
    )
  }

  out <- list(
    estimates = estimate,
    boundary = boundary,
    exposure = exposure,
    estimator = estimator,
    cutoff = cutoff,
    p = p,
    q = q,
    call = match.call()
  )

  class(out) <- "rddinterf"
  out
}


#' Summarize rddinterf estimation results
#'
#' Produces a compact summary table for objects returned by \code{\link{rddinterf}}.
#'
#' @param object An object of class \code{"rddinterf"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A numeric matrix with columns:
#'
#'   \describe{
#'     \item{`Estimate`}{Point estimate.}
#'     \item{`S.E.`}{Estimated standard error.}
#'     \item{`CI.l`}{Lower endpoint of the 95% confidence interval.}
#'     \item{`CI.u`}{Upper endpoint of the 95% confidence interval.}
#'     \item{`p.value`}{Two-sided p-value based on a normal approximation.}
#'     \item{`N.l`}{Effective sample size to the left of the cutoff.}
#'     \item{`N.r`}{Effective sample size to the right of the cutoff.}
#'   }
#'
#' @export
summary.rddinterf <- function(object, ...) {
  results <- object$estimates

  keep <- vapply(results, function(x) is.list(x) && !is.null(x$est), logical(1))
  results <- results[keep]

  if (length(results) == 0L) {
    out <- matrix(
      numeric(0),
      nrow = 0,
      ncol = 7
    )
    colnames(out) <- c("Estimate", "S.E.", "CI.l", "CI.u", "p.value", "N.l", "N.r")
    return(out)
  }

  out <- do.call(
    rbind,
    lapply(results, function(x) {
      c(
        Estimate = x$est,
        `S.E.` = x$se,
        `CI.l` = x$ci[1],
        `CI.u` = x$ci[2],
        `p.value` = 2 * (1 - pnorm(abs(x$est / x$se))),
        `N.l` = x$N_l,
        `N.r` = x$N_r
      )
    })
  )

  out <- as.matrix(out)
  rownames(out) <- names(results)
  out
}




#' Local polynomial estimator without bias correction
#'
#' Internal helper used by `rddinterf()` to estimate a cutoff contrast using a
#' local polynomial regression of order `p` and the network-robust variance
#' estimator under a supplied dependency graph.
#'
#' @param h Scalar estimation bandwidth.
#' @param X Numeric running variable.
#' @param Y Numeric outcome vector.
#' @param GG Dependency graph matrix for the sample.
#' @param subset Optional logical vector specifying the estimation sample.
#' @param p Integer polynomial order for the main local fit.
#'
#' @return A list containing the point estimate, estimated variance, effective
#'   sample sizes, side-specific variance components, fitted coefficient vectors,
#'   and side-specific subgraphs.
#'
#' @keywords internal
#' @noRd
estimate_conv <- function(h, X, Y, GG, subset = NULL, p = 1) {

  if (!is.null(subset)) {
    Y <- Y[subset]
    X <- X[subset]
    GG <- GG[subset, subset, drop = FALSE]
  }

  nas <- is.na(X) | is.na(Y)
  X <- X[!nas]
  Y <- Y[!nas]
  GG <- GG[!nas, !nas, drop = FALSE]

  n <- length(Y)
  if (n == 0L) {
    stop("No non-missing observations available.")
  }

  Z <- X >= 0

  # Polynomial bases:
  # Xp uses scaled running variable X / h
  # XX uses unscaled running variable X
  Xp <- outer(X / h, 0:p, `^`)
  XX <- outer(X, 0:p, `^`)

  # Triangular kernel weights
  omega_h_r <- Z * kernel1(X / h, h)
  omega_h_l <- (!Z) * kernel1(X / h, h)

  w_h_r <- omega_h_r > 0
  w_h_l <- omega_h_l > 0

  N_r <- sum(w_h_r)
  N_l <- sum(w_h_l)

  if (N_r == 0L || N_l == 0L) {
    stop("No observations on one side of the cutoff within bandwidth.")
  }

  GG_h_r <- GG[w_h_r, w_h_r, drop = FALSE]
  GG_h_l <- GG[w_h_l, w_h_l, drop = FALSE]
  GG_h_l_r <- GG[w_h_r, w_h_l, drop = FALSE]

  # Moment matrices
  gamma_h_r <- crossprod(Xp, Xp * omega_h_r) / n
  gamma_h_l <- crossprod(Xp, Xp * omega_h_l) / n

  gamma_h_r_inv <- solve(gamma_h_r)
  gamma_h_l_inv <- solve(gamma_h_l)

  # Rescaling matrix: diag(1, h^{-1}, ..., h^{-p})
  Hp <- diag((1 / h)^(0:p))

  # Coefficients on each side
  beta_plus <- Hp %*% gamma_h_r_inv %*% crossprod(Xp, omega_h_r * Y) / n
  beta_minus <- Hp %*% gamma_h_l_inv %*% crossprod(Xp, omega_h_l * Y) / n

  # Residuals
  res_r <- as.vector(Y - XX %*% beta_plus)
  res_l <- as.vector(Y - XX %*% beta_minus)

  # Score contributions
  M_h_r <- (omega_h_r * Xp * res_r)[w_h_r, , drop = FALSE]
  M_h_l <- (omega_h_l * Xp * res_l)[w_h_l, , drop = FALSE]

  psi_h_r <- crossprod(M_h_r, as.matrix(GG_h_r %*% M_h_r)) / n
  psi_h_l <- crossprod(M_h_l, as.matrix(GG_h_l %*% M_h_l)) / n
  psi_h_l_r <- crossprod(M_h_r, as.matrix(GG_h_l_r %*% M_h_l)) / n

  V_h_r <- Hp %*% gamma_h_r_inv %*% psi_h_r %*% gamma_h_r_inv %*% Hp / n
  V_h_l <- Hp %*% gamma_h_l_inv %*% psi_h_l %*% gamma_h_l_inv %*% Hp / n
  V_h_l_r <- Hp %*% gamma_h_r_inv %*% psi_h_l_r %*% gamma_h_l_inv %*% Hp / n

  # Intercept contrast at the cutoff
  V_r <- V_h_r[1, 1]
  V_l <- V_h_l[1, 1]
  C_r_l <- V_h_l_r[1, 1]
  V <- V_r + V_l - 2 * C_r_l

  tau <- beta_plus[1] - beta_minus[1]

  list(
    est = tau,
    var = V,
    N_r = N_r,
    N_l = N_l,
    V_r = V_r,
    V_l = V_l,
    C_r_l = C_r_l,
    beta_plus = beta_plus,
    beta_minus = beta_minus,
    GG_h_r = GG_h_r,
    GG_h_l = GG_h_l
  )
}

#' Bias-corrected local polynomial estimator
#'
#' Internal helper used by `rddinterf()` to estimate a cutoff contrast using a
#' local polynomial regression of order `p` with bias correction based on a
#' polynomial of order `q`, together with a network-robust variance estimator.
#'
#' @param h Scalar estimation bandwidth.
#' @param b Scalar bias-correction bandwidth.
#' @param X Numeric running variable.
#' @param Y Numeric outcome vector.
#' @param GG Dependency graph matrix for the sample.
#' @param subset Optional logical vector specifying the estimation sample.
#' @param p Integer polynomial order for the main local fit.
#' @param q Integer polynomial order for the bias-correction fit. Must satisfy
#'   `q >= p + 1`.
#'
#' @return A list containing the bias-corrected point estimate, estimated
#'   variance, effective sample sizes, and fitted coefficient vectors for the
#'   main and bias-correction regressions.
#'
#' @keywords internal
#' @noRd
estimate_bc <- function(h, b, X, Y, GG, subset = NULL, p = 1, q = p + 1) {
  if (q < p + 1) {
    stop("q must be at least p + 1.")
  }

  if (!is.null(subset)) {
    Y <- Y[subset]
    X <- X[subset]
    GG <- GG[subset, subset, drop = FALSE]
  }

  nas <- is.na(X) | is.na(Y)
  X <- X[!nas]
  Y <- Y[!nas]
  GG <- GG[!nas, !nas, drop = FALSE]

  n <- length(Y)
  if (n == 0L) {
    stop("No non-missing observations available.")
  }

  Z <- X >= 0

  # Polynomial bases
  Xp <- outer(X / h, 0:p, `^`)
  Xq <- outer(X / b, 0:q, `^`)
  XXq <- outer(X, 0:q, `^`)

  # Kernel weights
  omega_h_r <- Z * kernel1(X / h, h)
  omega_h_l <- (!Z) * kernel1(X / h, h)

  omega_b_r <- Z * kernel1(X / b, b)
  omega_b_l <- (!Z) * kernel1(X / b, b)

  w_h_r <- omega_h_r > 0
  w_h_l <- omega_h_l > 0
  w_b_r <- omega_b_r > 0
  w_b_l <- omega_b_l > 0

  N_r <- sum(w_h_r)
  N_l <- sum(w_h_l)

  if (N_r == 0L || N_l == 0L) {
    stop("No observations on one side of the cutoff within estimation bandwidth.")
  }

  # Restrict objects by side/bandwidth
  Xp_h_r <- Xp[w_h_r, , drop = FALSE]
  Xp_h_l <- Xp[w_h_l, , drop = FALSE]
  Xq_b_r <- Xq[w_b_r, , drop = FALSE]
  Xq_b_l <- Xq[w_b_l, , drop = FALSE]

  omega_h_r_sub <- omega_h_r[w_h_r]
  omega_h_l_sub <- omega_h_l[w_h_l]
  omega_b_r_sub <- omega_b_r[w_b_r]
  omega_b_l_sub <- omega_b_l[w_b_l]

  Y_h_r <- Y[w_h_r]
  Y_h_l <- Y[w_h_l]
  Y_b_r <- Y[w_b_r]
  Y_b_l <- Y[w_b_l]

  GG_h_r <- GG[w_h_r, w_h_r, drop = FALSE]
  GG_h_l <- GG[w_h_l, w_h_l, drop = FALSE]
  GG_h_l_r <- GG[w_h_r, w_h_l, drop = FALSE]

  GG_b_r <- GG[w_b_r, w_b_r, drop = FALSE]
  GG_b_l <- GG[w_b_l, w_b_l, drop = FALSE]
  GG_b_l_r <- GG[w_b_r, w_b_l, drop = FALSE]

  GG_h_b_r_r <- GG[w_h_r, w_b_r, drop = FALSE]
  GG_h_b_l_l <- GG[w_h_l, w_b_l, drop = FALSE]
  GG_h_b_r_l <- GG[w_h_r, w_b_l, drop = FALSE]
  GG_h_b_l_r <- GG[w_b_r, w_h_l, drop = FALSE]

  # Scaling matrices
  Hp <- diag((1 / h)^(0:p))
  Hq <- diag((1 / b)^(0:q))

  # Moment matrices
  gamma_h_r <- crossprod(Xp_h_r, Xp_h_r * omega_h_r_sub) / n
  gamma_h_l <- crossprod(Xp_h_l, Xp_h_l * omega_h_l_sub) / n
  gamma_b_r <- crossprod(Xq_b_r, Xq_b_r * omega_b_r_sub) / n
  gamma_b_l <- crossprod(Xq_b_l, Xq_b_l * omega_b_l_sub) / n

  gamma_h_r_inv <- solve(gamma_h_r)
  gamma_h_l_inv <- solve(gamma_h_l)
  gamma_b_r_inv <- solve(gamma_b_r)
  gamma_b_l_inv <- solve(gamma_b_l)

  # Main p-order estimators
  beta_plus <- Hp %*% gamma_h_r_inv %*% crossprod(Xp_h_r, omega_h_r_sub * Y_h_r) / n
  beta_minus <- Hp %*% gamma_h_l_inv %*% crossprod(Xp_h_l, omega_h_l_sub * Y_h_l) / n

  # q-order estimators used for bias correction
  beta_plus_q <- Hq %*% gamma_b_r_inv %*% crossprod(Xq_b_r, omega_b_r_sub * Y_b_r) / n
  beta_minus_q <- Hq %*% gamma_b_l_inv %*% crossprod(Xq_b_l, omega_b_l_sub * Y_b_l) / n

  # Bias constants
  # theta = E[ r_p(u) K(u) u^(p+1) ]
  theta_p_q_r <- crossprod(Xp_h_r, omega_h_r_sub * (X[w_h_r] / h)^(p + 1)) / n
  theta_p_q_l <- crossprod(Xp_h_l, omega_h_l_sub * (X[w_h_l] / h)^(p + 1)) / n

  B_r <- gamma_h_r_inv %*% theta_p_q_r
  B_l <- gamma_h_l_inv %*% theta_p_q_l

  # Residuals from q-order fit
  res_r <- as.vector(Y - XXq %*% beta_plus_q)
  res_l <- as.vector(Y - XXq %*% beta_minus_q)

  # Score contributions for h-bandwidth part
  M_h_r <- (omega_h_r_sub * Xp_h_r) * res_r[w_h_r]
  M_h_l <- (omega_h_l_sub * Xp_h_l) * res_l[w_h_l]

  psi_h_r <- crossprod(M_h_r, GG_h_r %*% M_h_r) / n
  psi_h_l <- crossprod(M_h_l, GG_h_l %*% M_h_l) / n
  psi_h_l_r <- crossprod(M_h_r, GG_h_l_r %*% M_h_l) / n

  V_h_r <- Hp %*% gamma_h_r_inv %*% psi_h_r %*% gamma_h_r_inv %*% Hp / n
  V_h_l <- Hp %*% gamma_h_l_inv %*% psi_h_l %*% gamma_h_l_inv %*% Hp / n
  V_h_l_r <- Hp %*% gamma_h_r_inv %*% psi_h_l_r %*% gamma_h_l_inv %*% Hp / n

  # Score contributions for b-bandwidth part
  M_b_r <- (omega_b_r_sub * Xq_b_r) * res_r[w_b_r]
  M_b_l <- (omega_b_l_sub * Xq_b_l) * res_l[w_b_l]

  psi_b_r <- crossprod(M_b_r, GG_b_r %*% M_b_r) / n
  psi_b_l <- crossprod(M_b_l, GG_b_l %*% M_b_l) / n
  psi_b_l_r <- crossprod(M_b_r, GG_b_l_r %*% M_b_l) / n

  V_b_r <- Hq %*% gamma_b_r_inv %*% psi_b_r %*% gamma_b_r_inv %*% Hq / n
  V_b_l <- Hq %*% gamma_b_l_inv %*% psi_b_l %*% gamma_b_l_inv %*% Hq / n
  V_b_l_r <- Hq %*% gamma_b_r_inv %*% psi_b_l_r %*% gamma_b_l_inv %*% Hq / n

  # Cross-covariances between h and b parts
  psi_p_q_r <- crossprod(M_h_r, GG_h_b_r_r %*% M_b_r) / n
  psi_p_q_l <- crossprod(M_h_l, GG_h_b_l_l %*% M_b_l) / n

  C_p_q_r <- Hp %*% gamma_h_r_inv %*% psi_p_q_r %*% gamma_b_r_inv %*% Hq / n
  C_p_q_l <- Hp %*% gamma_h_l_inv %*% psi_p_q_l %*% gamma_b_l_inv %*% Hq / n

  psi_p_q_l_r <- crossprod(M_h_r, GG_h_b_r_l %*% M_b_l) / n
  psi_p_q_l_r2 <- crossprod(M_b_r, GG_h_b_l_r %*% M_h_l) / n

  C_p_q_l_r <- Hp %*% gamma_h_r_inv %*% psi_p_q_l_r %*% gamma_b_l_inv %*% Hq / n
  C_p_q_l_r2 <- Hq %*% gamma_b_r_inv %*% psi_p_q_l_r2 %*% gamma_h_l_inv %*% Hp / n

  # Index of x^(p+1) coefficient in q-order fit
  idx_bias <- p + 2

  # Robust variance pieces
  V_r_bc <- V_h_r[1, 1] +
    h^(2 * (p + 1)) * V_b_r[idx_bias, idx_bias] * (B_r[1, 1]^2) -
    2 * h^(p + 1) * C_p_q_r[1, idx_bias] * B_r[1, 1]

  V_l_bc <- V_h_l[1, 1] +
    h^(2 * (p + 1)) * V_b_l[idx_bias, idx_bias] * (B_l[1, 1]^2) -
    2 * h^(p + 1) * C_p_q_l[1, idx_bias] * B_l[1, 1]

  C_bc <- V_h_l_r[1, 1] +
    h^(2 * (p + 1)) * V_b_l_r[idx_bias, idx_bias] * B_r[1, 1] * B_l[1, 1] -
    h^(p + 1) * C_p_q_l_r[1, idx_bias] * B_r[1, 1] -
    h^(p + 1) * C_p_q_l_r2[idx_bias, 1] * B_l[1, 1]

  V_bc <- V_r_bc + V_l_bc - 2 * C_bc

  # Bias-corrected side-specific limits and their contrast
  beta_plus_bc <-
    beta_plus[1] - h^(p + 1) * beta_plus_q[idx_bias] * B_r[1, 1]
  beta_minus_bc <-
    beta_minus[1] - h^(p + 1) * beta_minus_q[idx_bias] * B_l[1, 1]
  tau_bc <- beta_plus_bc - beta_minus_bc

  list(
    est = tau_bc,
    var = V_bc,
    N_r = N_r,
    N_l = N_l,
    beta_plus = beta_plus,
    beta_minus = beta_minus,
    beta_plus_q = beta_plus_q,
    beta_minus_q = beta_minus_q,
    beta_plus_bc = beta_plus_bc,
    beta_minus_bc = beta_minus_bc
  )
}





#' Fit a single estimand
#'
#' Internal helper used by `rddinterf()` to estimate one overall or
#' contrast-specific effect, select bandwidths, and return a standardized result
#' object.
#'
#' @param name Character name of the estimand.
#' @param Xrun Numeric running variable used for estimation.
#' @param keep Logical vector indicating the estimation sample.
#' @param boundary_obj Boundary object or boundary description associated with
#'   the estimand.
#' @param sbar Boundary codimension used by the bandwidth selector.
#' @param Y Numeric outcome vector.
#' @param GG Dependency graph matrix for the sample.
#' @param estimator Character string specifying `"conventional"` or `"robust"`.
#' @param p Integer polynomial order for the main local fit.
#' @param q Integer polynomial order for the bias-correction fit.
#' @param h_manual Optional manually supplied estimation bandwidth.
#' @param b_manual Optional manually supplied bias-correction bandwidth.
#'
#' @return Either `NA` if estimation fails, or a list containing the point
#'   estimate, standard error, confidence interval, boundary object, effective
#'   sample sizes, bandwidths, estimation sample indicator, and polynomial
#'   orders.
#'
#' @keywords internal
#' @noRd
fit_one <- function(
    name,
    Xrun,
    keep,
    boundary_obj,
    sbar,
    Y,
    GG,
    estimator,
    p = 1,
    q = p + 1,
    h_manual = NULL,
    b_manual = NULL)
{
  if (q < p + 1) {
    stop("q must be at least p + 1.")
  }

  bws_auto <- get_auto_bw(
    Y = Y,
    Xrun = Xrun,
    GG = GG,
    subset = keep,
    sbar = sbar,
    p = p,
    q = q
  )

  bw_use <- resolve_bw(name, h_manual, b_manual, bws_auto)

  supplied <- !is.na(bw_use)
  if (any(supplied & (!is.finite(bw_use) | bw_use <= 0))) {
    stop("Bandwidths must be positive finite numbers.")
  }

  # -----------------------------
  # Conventional estimator
  # -----------------------------
  if (estimator == "conventional") {
    if (is.na(bw_use["h"])) return(NA)

    fit <- tryCatch(
      estimate_conv(
        h = bw_use["h"],
        X = Xrun,
        Y = Y,
        GG = GG,
        subset = keep,
        p = p
      ),
      error = function(e) {
        message("estimate_conv failed for ", name, ": ", e$message)
        NULL
      }
    )

    if (is.null(fit)) return(NA)

    lim <- 1.96 * sqrt(fit$var)

    return(list(
      est = fit$est,
      se = sqrt(fit$var),
      ci = fit$est + c(-lim, lim),
      mean_left = as.numeric(fit$beta_minus[1]),
      mean_right = as.numeric(fit$beta_plus[1]),
      boundary = boundary_obj,
      dstar = Xrun[keep],
      Yi = Y[keep],
      N_r = fit$N_r,
      N_l = fit$N_l,
      h = bw_use["h"],
      b = bw_use["b"],
      tr = keep,
      p = p,
      q = q
    ))
  }

  # -----------------------------
  # Robust estimator
  # -----------------------------
  if (is.na(bw_use["h"]) || is.na(bw_use["b"])) return(NA)

  fit <- tryCatch(
    estimate_bc(
      h = bw_use["h"],
      b = bw_use["b"],
      X = Xrun,
      Y = Y,
      GG = GG,
      subset = keep,
      p = p,
      q = q
    ),
    error = function(e) {
      message("estimate_bc failed for ", name, ": ", e$message)
      NULL
    }
  )


  if (is.null(fit)) return(NA)

  lim <- 1.96 * sqrt(fit$var)

  list(
    est = fit$est,
    se = sqrt(fit$var),
    ci = fit$est + c(-lim, lim),
    mean_left = as.numeric(fit$beta_minus_bc),
    mean_right = as.numeric(fit$beta_plus_bc),
    boundary = boundary_obj,
    dstar = Xrun[keep],
    Yi = Y[keep],
    N_r = fit$N_r,
    N_l = fit$N_l,
    h = bw_use["h"],
    b = bw_use["b"],
    tr = keep,
    p = p,
    q = q
  )
}


# kernel1: triangular kernel function.
kernel1 = function(x, h){

  k = 1 - abs(x)
  k[k < 0 ] = 0
  k = k/h
  return(k)

}











