#' Distance to the boundary under one-treated exposure
#'
#' Internal helper that computes the Euclidean distance from each observation to
#' the separating boundary between two effective treatments under the
#' `"one.treated"` exposure mapping.
#'
#' @param Xi Numeric vector of centered own running variables.
#' @param Xj List of numeric neighbor running-variable vectors.
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param no_neighbor_value Value returned for units with no neighbors when the
#'   exposure contrast requires a change in neighbor treatment status.
#'
#' @return Numeric vector of distances, with one entry per observation.
#'
#' @keywords internal
#' @noRd
min_distance_one_treated <- function(
    Xi, Xj, d, g, dprime, gprime,
    no_neighbor_value = NA_real_
) {
  stopifnot(is.numeric(Xi), is.list(Xj), length(Xi) == length(Xj))

  n <- length(Xi)
  mlen <- vapply(Xj, length, integer(1))
  dims <- sort(unique(mlen))

  out <- rep(NA_real_, n)

  need_x <- (d != dprime)
  need_g <- (g != gprime)

  for (m in dims) {
    idx <- which(mlen == m)
    if (!length(idx)) next

    dx2 <- if (need_x) Xi[idx]^2 else rep(0, length(idx))

    if (!need_g) {
      out[idx] <- sqrt(dx2)
      next
    }

    if (m == 0L) {
      out[idx] <- no_neighbor_value
      next
    }

    V <- matrix(unlist(Xj[idx]), ncol = m, byrow = TRUE)

    M <- apply(V, 1, max)
    vpos2 <- rowSums((V * (V > 0))^2)

    dg2 <- ifelse(M > 0, vpos2, M^2)

    out[idx] <- sqrt(dx2 + dg2)
  }

  out
}

#' Distance to the boundary under sum-treated exposure
#'
#' Internal helper that computes the Euclidean distance from each observation to
#' the separating boundary between two effective treatments under the
#' `"sum.treated"` exposure mapping.
#'
#' @param Xi Numeric vector of centered own running variables.
#' @param Xj List of numeric neighbor running-variable vectors.
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param no_neighbor_value Value returned for units with no neighbors when the
#'   exposure contrast requires a change in neighbor treatment status.
#'
#' @return Numeric vector of distances, with one entry per observation.
#'
#' @keywords internal
#' @noRd
min_distance_sum_treated <- function(
    Xi, Xj, d, g, dprime, gprime,
    no_neighbor_value = NA_real_
) {
  stopifnot(is.numeric(Xi), is.list(Xj), length(Xi) == length(Xj))

  n <- length(Xi)
  mlen <- vapply(Xj, length, integer(1))
  dims <- sort(unique(mlen))

  need_x <- (d != dprime)
  need_g <- (g != gprime)

  out <- rep(NA_real_, n)

  for (m in dims) {
    idx <- which(mlen == m)
    if (!length(idx)) next

    dx2 <- if (need_x) Xi[idx]^2 else rep(0, length(idx))

    if (!need_g) {
      out[idx] <- sqrt(dx2)
      next
    }

    if (m == 0L) {
      out[idx] <- no_neighbor_value
      next
    }

    dg <- vapply(
      idx,
      function(k) dist_sum_treated_boundary_one(Xj[[k]], g, gprime),
      numeric(1)
    )

    out[idx] <- sqrt(dx2 + dg^2)
  }

  out
}

#' Distance to the boundary under fraction-treated exposure
#'
#' Internal helper that computes the Euclidean distance from each observation to
#' the separating boundary between two effective treatments under the
#' `"frac.treated"` exposure mapping.
#'
#' @param Xi Numeric vector of centered own running variables.
#' @param Xj List of numeric neighbor running-variable vectors.
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param no_neighbor_value Value returned for units with no neighbors when the
#'   exposure contrast requires a change in neighbor treatment status.
#' @param tol Tolerance used to determine whether `g * m` and `gprime * m` are
#'   integer-compatible for neighborhood size `m`.
#'
#' @return Numeric vector of distances, with one entry per observation.
#'
#' @keywords internal
#' @noRd
min_distance_frac_treated <- function(
    Xi, Xj, d, g, dprime, gprime,
    no_neighbor_value = NA_real_,
    tol = 1e-8
) {
  stopifnot(is.numeric(Xi), is.list(Xj), length(Xi) == length(Xj))

  n <- length(Xi)
  mlen <- vapply(Xj, length, integer(1))
  dims <- sort(unique(mlen))

  need_x <- (d != dprime)
  need_g <- (g != gprime)

  out <- rep(NA_real_, n)

  for (m in dims) {
    idx <- which(mlen == m)
    if (!length(idx)) next

    dx2 <- if (need_x) Xi[idx]^2 else rep(0, length(idx))

    if (!need_g) {
      out[idx] <- sqrt(dx2)
      next
    }

    if (m == 0L) {
      out[idx] <- no_neighbor_value
      next
    }

    k <- g * m
    kp <- gprime * m

    if (abs(k - round(k)) > tol || abs(kp - round(kp)) > tol) {
      out[idx] <- NA_real_
      next
    }

    k <- as.integer(round(k))
    kp <- as.integer(round(kp))

    dg <- vapply(
      idx,
      function(j) dist_sum_treated_boundary_one(Xj[[j]], k, kp),
      numeric(1)
    )

    out[idx] <- sqrt(dx2 + dg^2)
  }

  out
}

#' Distance to the boundary under majority exposure
#'
#' Internal helper that computes the Euclidean distance from each observation to
#' the separating boundary between two effective treatments under the
#' `"majority"` exposure mapping.
#'
#' @param Xi Numeric vector of centered own running variables.
#' @param Xj List of numeric neighbor running-variable vectors.
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param no_neighbor_value Value returned for units with no neighbors when the
#'   exposure contrast requires a change in neighbor treatment status.
#'
#' @return Numeric vector of distances, with one entry per observation.
#'
#' @keywords internal
#' @noRd
min_distance_majority <- function(
    Xi, Xj, d, g, dprime, gprime,
    no_neighbor_value = NA_real_
) {
  stopifnot(is.numeric(Xi), is.list(Xj), length(Xi) == length(Xj))

  n <- length(Xi)
  mlen <- vapply(Xj, length, integer(1))
  dims <- sort(unique(mlen))

  need_x <- (d != dprime)
  need_g <- (g != gprime)

  out <- rep(NA_real_, n)

  for (m in dims) {
    idx <- which(mlen == m)
    if (!length(idx)) next

    dx2 <- if (need_x) Xi[idx]^2 else rep(0, length(idx))

    if (!need_g) {
      out[idx] <- sqrt(dx2)
      next
    }

    if (m == 0L) {
      out[idx] <- no_neighbor_value
      next
    }

    q <- ceiling(m / 2)

    dg <- vapply(
      idx,
      function(j) dist_sum_treated_boundary_one(Xj[[j]], q, q - 1L),
      numeric(1)
    )

    out[idx] <- sqrt(dx2 + dg^2)
  }

  out
}

#' Dispatch minimum-distance computation by exposure mapping
#'
#' Internal helper used by `compute_dstar()` to select the appropriate
#' minimum-distance routine for the requested exposure mapping.
#'
#' @param exposure Character string specifying the exposure mapping.
#' @param Xi Numeric vector of centered own running variables.
#' @param Xj List of numeric neighbor running-variable vectors.
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#'
#' @return Numeric vector of distances for the requested contrast.
#'
#' @keywords internal
#' @noRd
compute_min_distance <- function(exposure, Xi, Xj, d, g, dprime, gprime) {
  if (exposure == "one.treated") {
    return(min_distance_one_treated(Xi, Xj, d, g, dprime, gprime))
  }

  if (exposure == "sum.treated") {
    return(min_distance_sum_treated(Xi, Xj, d, g, dprime, gprime))
  }

  if (exposure == "frac.treated") {
    return(min_distance_frac_treated(Xi, Xj, d, g, dprime, gprime))
  }

  if (exposure == "majority") {
    return(min_distance_majority(Xi, Xj, d, g, dprime, gprime))
  }

  stop("Unsupported exposure.")
}

#' Construct the signed running variable for a treatment contrast
#'
#' Internal helper that computes the distance-based running variable used for a
#' contrast-specific estimand under the selected exposure mapping.
#'
#' @param Xi Numeric vector of centered own running variables.
#' @param Xj List of numeric neighbor running-variable vectors.
#' @param keep Logical vector indicating the observations used for the contrast.
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param exposure Character string specifying the exposure mapping.
#'
#' @return Numeric vector of signed running variables, with `NA` outside the
#'   estimation sample.
#'
#' @keywords internal
#' @noRd
compute_dstar <- function(
    Xi, Xj, keep,
    d, g, dprime, gprime,
    exposure
) {
  dstar <- rep(NA_real_, length(Xi))

  if (g == gprime) {
    dstar[keep] <- Xi[keep]
    return(dstar)
  }

  dtmp <- compute_min_distance(
    exposure = exposure,
    Xi = Xi[keep],
    Xj = Xj[keep],
    d = d,
    g = g,
    dprime = dprime,
    gprime = gprime
  )

  dstar[keep] <- dtmp

  dstar[keep & (dstar < 0)] <- -abs(dstar[keep & (dstar < 0)])
  dstar[keep & (dstar > 0)] <-  abs(dstar[keep & (dstar > 0)])

  dstar
}

#' Distance to a sum-treated separating boundary for one neighbor-score vector
#'
#' Internal helper used by the distance routines for `"sum.treated"`,
#' `"frac.treated"`, and `"majority"` exposures. For one vector of neighbor
#' running variables, it computes the Euclidean distance to the closest
#' separating boundary between two treated-neighbor counts.
#'
#' @param v Numeric vector of neighbor running variables.
#' @param g Integer treated-neighbor count for one effective treatment state.
#' @param gprime Integer treated-neighbor count for the comparison treatment state.
#'
#' @return A nonnegative numeric scalar giving the Euclidean distance from `v`
#'   to the nearest separating boundary, or `NA_real_` if the input is not
#'   compatible with the requested counts.
#'
#' @keywords internal
#' @noRd
dist_sum_treated_boundary_one <- function(v, g, gprime) {
  m <- length(v)
  if (m == 0L) return(NA_real_)

  g <- as.integer(g)
  gprime <- as.integer(gprime)

  if (g < 0L || gprime < 0L || g > m || gprime > m) {
    return(NA_real_)
  }

  gl <- min(g, gprime)
  gh <- max(g, gprime)
  r <- gh - gl

  if (r == 0L) return(0)

  a <- sort(v, decreasing = TRUE)

  cost <- 0

  cost <- cost + sum(a[(gl + 1L):gh]^2)

  if (gl > 0L) {
    ap <- a[1L:gl]
    cost <- cost + sum(ap[ap < 0]^2)
  }

  if (gh < m) {
    an <- a[(gh + 1L):m]
    cost <- cost + sum(an[an > 0]^2)
  }

  sqrt(cost)
}
