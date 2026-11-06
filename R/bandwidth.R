#' Network-robust bandwidth selector
#'
#' Internal helper that computes automatic estimation and bias-correction
#' bandwidths for the local polynomial estimator under interference.
#'
#' The selector uses pilot quantities based on `get_var_bias()` and returns the
#' bandwidth vector in the same four-entry format used throughout the package:
#' left and right estimation bandwidths followed by left and right
#' bias-correction bandwidths.
#'
#' @param X1 Numeric running variable.
#' @param Y1 Numeric outcome vector.
#' @param p Integer polynomial order for the main local fit.
#' @param q Integer polynomial order for the bias-correction fit. Must satisfy
#'   `q >= p + 1`.
#' @param sbar Boundary codimension used in the bandwidth formula.
#' @param GG Dependency graph matrix for the sample.
#' @param subset Optional logical vector specifying the estimation sample.
#' @param masspoints Logical; if `TRUE`, use the number of unique running-variable
#'   values in the pilot bandwidth calculation.
#'
#' @return Numeric vector of four bandwidths:
#'   `c(h_left, h_right, b_left, b_right)`.
#'
#' @keywords internal
#' @noRd
bws_robust <- function(X1, Y1, p = 1, q = p + 1, sbar, GG,
                       subset = NULL, masspoints = TRUE) {
  if (q < p + 1) {
    stop("q must be at least p + 1.")
  }

  if (!is.null(subset)) {
    X1 <- X1[subset]
    Y1 <- Y1[subset]
    GG <- GG[subset, subset, drop = FALSE]
  }

  nas <- is.na(Y1) | is.na(X1)
  X1 <- X1[!nas]
  Y1 <- Y1[!nas]
  GG <- GG[!nas, !nas, drop = FALSE]

  Z1 <- X1 >= 0
  n <- length(X1)

  if (n == 0L) {
    stop("No non-missing observations available.")
  }

  if (!any(Z1) || !any(!Z1)) {
    stop("Need observations on both sides of the cutoff.")
  }

  scale_x <- min(
    sd(X1),
    (quantile(X1, 0.75, type = 2) - quantile(X1, 0.25, type = 2)) / 1.349
  )

  if (masspoints) {
    M <- length(unique(X1))
    vn <- 2.576 * scale_x * M^(-1 / (4 + sbar))
  } else {
    vn <- 2.576 * scale_x * n^(-1 / (4 + sbar))
  }

  r_l <- abs(min(X1[!Z1]))
  r_r <- abs(max(X1[Z1]))

  X1_l <- X1[!Z1]
  Y1_l <- Y1[!Z1]
  GG_l <- GG[!Z1, !Z1, drop = FALSE]

  X1_r <- X1[Z1]
  Y1_r <- Y1[Z1]
  GG_r <- GG[Z1, Z1, drop = FALSE]

  bw_max <- max(r_l, r_r)

  C_d_l1 <- get_var_bias(
    X1_l, Y1_l,
    h_v = vn, h_B = r_l,
    GG = GG_l,
    o = q + 1, nu = q + 1, o_B = q + 2,
    ss = sbar
  )
  C_d_r1 <- get_var_bias(
    X1_r, Y1_r,
    h_v = vn, h_B = r_r,
    GG = GG_r,
    o = q + 1, nu = q + 1, o_B = q + 2,
    ss = sbar
  )

  d_bw_d1 <- ((C_d_l1$V + C_d_r1$V) / (C_d_r1$B - C_d_l1$B)^2)^C_d_l1$rate
  d_bw_d1 <- min(d_bw_d1, bw_max)

  C_b_l1 <- get_var_bias(
    X1_l, Y1_l,
    h_v = vn, h_B = d_bw_d1,
    GG = GG_l,
    o = q, nu = p + 1, o_B = q + 1,
    ss = sbar, scale = 1
  )
  C_b_r1 <- get_var_bias(
    X1_r, Y1_r,
    h_v = vn, h_B = d_bw_d1,
    GG = GG_r,
    o = q, nu = p + 1, o_B = q + 1,
    ss = sbar, scale = 1
  )

  b_bw_d1 <- ((C_b_l1$V + C_b_r1$V) /
                ((C_b_r1$B - C_b_l1$B)^2 + (C_b_r1$R + C_b_l1$R)))^C_b_l1$rate
  b_bw_d1 <- min(b_bw_d1, bw_max)

  C_h_l1 <- get_var_bias(
    X1_l, Y1_l,
    h_v = vn, h_B = b_bw_d1,
    GG = GG_l,
    o = p, nu = 0, o_B = q,
    ss = sbar, scale = 1
  )
  C_h_r1 <- get_var_bias(
    X1_r, Y1_r,
    h_v = vn, h_B = b_bw_d1,
    GG = GG_r,
    o = p, nu = 0, o_B = q,
    ss = sbar, scale = 1
  )

  h_bw_d1 <- ((C_h_l1$V + C_h_r1$V) /
                ((C_h_r1$B - C_h_l1$B)^2 + (C_h_r1$R + C_h_l1$R)))^C_h_l1$rate
  h_bw_d1 <- min(h_bw_d1, bw_max)

  c(h_bw_d1, h_bw_d1, b_bw_d1, b_bw_d1)
}

#' Bias and variance building blocks for bandwidth selection
#'
#' Internal helper used by `bws_robust()` to compute the leading variance,
#' bias, and regularization terms entering the network-robust bandwidth
#' formulas for a given local polynomial order and derivative order.
#'
#' @param X Numeric running variable.
#' @param Y Numeric outcome vector.
#' @param h_v Pilot bandwidth used for the variance component.
#' @param h_B Pilot bandwidth used for the bias component.
#' @param GG Dependency graph matrix for the sample.
#' @param o Integer polynomial order used in the variance component.
#' @param nu Integer derivative order of interest.
#' @param o_B Integer polynomial order used in the bias component.
#' @param ss Boundary codimension.
#' @param scale Numeric regularization scale used in the bandwidth formula.
#'
#' @return A list with components:
#' - `V`: leading variance term.
#' - `B`: leading bias term.
#' - `R`: regularization term.
#' - `rate`: rate exponent used in the bandwidth formula.
#'
#' @keywords internal
#' @noRd
get_var_bias <- function(X, Y, h_v, h_B, GG, o, nu, o_B, ss, scale = 0) {
  n <- length(X)

  if (n == 0L) {
    stop("No observations available in get_var_bias().")
  }

  max_o <- max(o, o_B)
  Xpow <- outer(X, 0:max_o, `^`)

  Xq <- sweep(Xpow[, 1:(o + 1), drop = FALSE], 2, h_v^(0:o), `/`)
  XX <- Xpow[, 1:(o + 1), drop = FALSE]
  w_v <- kernel1(X / h_v, h_v)

  gamma_hv <- crossprod(Xq, Xq * w_v) / n
  gamma_hv_inv <- solve(gamma_hv)

  Hq <- diag((1 / h_v)^(0:o))

  XtWY <- crossprod(Xq, w_v * Y) / n
  gamma_v <- Hq %*% (gamma_hv_inv %*% XtWY)

  res_gamma_v <- as.vector(Y - XX %*% gamma_v)

  M_v <- Xq * (w_v * res_gamma_v)
  psi_v <- crossprod(M_v, GG %*% M_v) / n

  V_v_mat <- Hq %*% gamma_hv_inv %*% psi_v %*% gamma_hv_inv %*% Hq / n
  V_v <- V_v_mat[nu + 1, nu + 1]

  Xq_B <- sweep(Xpow[, 1:(o_B + 1), drop = FALSE], 2, h_B^(0:o_B), `/`)
  XX_B <- Xpow[, 1:(o_B + 1), drop = FALSE]
  w_b <- kernel1(X / h_B, h_B)

  gamma_hb <- crossprod(Xq_B, Xq_B * w_b) / n
  gamma_hb_inv <- solve(gamma_hb)

  Hq_b <- diag((1 / h_B)^(0:o_B))

  XtWY_b <- crossprod(Xq_B, w_b * Y) / n
  gamma_b <- Hq_b %*% (gamma_hb_inv %*% XtWY_b)

  res_gamma_b <- as.vector(Y - XX_B %*% gamma_b)

  M_b <- Xq_B * (w_b * res_gamma_b)
  psi_b <- crossprod(M_b, GG %*% M_b) / n

  V_b_mat <- Hq_b %*% gamma_hb_inv %*% psi_b %*% gamma_hb_inv %*% Hq_b / n

  idx_bias <- o + 2
  V_b <- V_b_mat[idx_bias, idx_bias]

  theta_v <- crossprod(Xq, w_v * (X / h_v)^(o + 1)) / n
  B_v <- (gamma_hv_inv %*% theta_v)[nu + 1]

  BWreg <- 3 * B_v^2 * V_b

  B <- sqrt(2 * (o + 1 - nu)) * B_v * gamma_b[idx_bias]
  V <- (2 * nu + ss) * h_v^(2 * nu + ss) * V_v
  R <- scale * (2 * (o + 1 - nu)) * BWreg
  rate <- 1 / (2 * o + ss + 2)

  list(V = V, B = B, R = R, rate = rate)
}

#' Compute automatic bandwidths
#'
#' Internal helper that computes the automatic estimation and bias-correction
#' bandwidths used by `rddinterf()`.
#'
#' @param Y Numeric outcome vector.
#' @param Xrun Numeric running variable used for a given estimand.
#' @param GG Dependency graph matrix for the sample.
#' @param subset Optional logical vector specifying the estimation sample.
#' @param sbar Boundary codimension used by the bandwidth selector.
#' @param p Integer polynomial order for the main local fit.
#' @param q Integer polynomial order for the bias-correction fit.
#'
#' @return A numeric vector of bandwidths, or `NULL` if bandwidth selection fails.
#'
#' @keywords internal
#' @noRd
get_auto_bw <- function(Y, Xrun, GG, subset, sbar, p = 1, q = p + 1) {
  tryCatch(
    bws_robust(
      X1 = Xrun,
      Y1 = Y,
      p = p,
      q = q,
      sbar = sbar,
      GG = GG,
      subset = subset,
      masspoints = TRUE
    ),
    error = function(e) NULL
  )
}


#' Resolve bandwidth inputs for a given estimand
#'
#' Internal helper that determines the estimation (`h`) and bias-correction (`b`)
#' bandwidths for a given estimand, combining automatically computed bandwidths
#' with optional user-supplied overrides.
#'
#' Manual bandwidths can be provided as:
#' - a single numeric value (applied to all estimands),
#' - a named numeric vector, or
#' - a named list with entries corresponding to estimand names.
#'
#' @param name Character name of the estimand.
#' @param h_manual Optional user-supplied estimation bandwidth(s).
#' @param b_manual Optional user-supplied bias-correction bandwidth(s).
#' @param bws_auto Optional numeric vector of automatic bandwidths, assumed to
#'   be of the form `c(h_left, h_right, b_left, b_right)`.
#'
#' @return Named numeric vector `c(h = ..., b = ...)`.
#'
#' @keywords internal
#' @noRd
resolve_bw <- function(name, h_manual = NULL, b_manual = NULL, bws_auto = NULL) {

  # Default: use automatic bandwidths if available
  h_out <- if (!is.null(bws_auto)) bws_auto[1] else NA_real_
  b_out <- if (!is.null(bws_auto)) bws_auto[3] else NA_real_

  extract_manual <- function(x, name) {
    if (is.null(x)) return(NULL)

    # scalar → apply to all
    if (is.numeric(x) && length(x) == 1 && is.null(names(x))) {
      return(x)
    }

    # named list
    if (is.list(x) && !is.null(x[[name]])) {
      return(x[[name]])
    }

    # named vector
    if (is.numeric(x) && !is.null(names(x)) && name %in% names(x)) {
      return(x[[name]])
    }

    NULL
  }

  h_val <- extract_manual(h_manual, name)
  b_val <- extract_manual(b_manual, name)

  if (!is.null(h_val)) h_out <- h_val
  if (!is.null(b_val)) b_out <- b_val

  c(h = unname(h_out), b = unname(b_out))
}

