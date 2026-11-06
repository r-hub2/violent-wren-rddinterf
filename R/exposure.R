#' Neighborhood treatment under one-treated exposure
#'
#' Internal helper that computes the neighborhood treatment indicator under the
#' `"one.treated"` exposure mapping.
#'
#' @param Di Numeric or integer vector of own treatment assignments.
#' @param neighbors List of neighbor indices for each observation.
#'
#' @return Numeric vector of neighborhood treatment values.
#'
#' @keywords internal
#' @noRd
one_treated <- function(Di, neighbors) {
  as.numeric(sapply(neighbors, function(j) !all(Di[j] == 0)))
}

#' Neighborhood treatment under sum-treated exposure
#'
#' Internal helper that computes the number of treated neighbors for each
#' observation.
#'
#' @param Di Numeric or integer vector of own treatment assignments.
#' @param neighbors List of neighbor indices for each observation.
#'
#' @return Numeric vector of treated-neighbor counts.
#'
#' @keywords internal
#' @noRd
sum_treated <- function(Di, neighbors) {
  as.numeric(sapply(neighbors, function(j) sum(Di[j])))
}

#' Neighborhood treatment under fraction-treated exposure
#'
#' Internal helper that computes the fraction of treated neighbors for each
#' observation.
#'
#' @param Di Numeric or integer vector of own treatment assignments.
#' @param neighbors List of neighbor indices for each observation.
#'
#' @return Numeric vector of treated-neighbor fractions.
#'
#' @keywords internal
#' @noRd
frac_treated <- function(Di, neighbors) {
  as.numeric(sapply(neighbors, function(j) {
    if (length(j) == 0L) return(NA_real_)
    sum(Di[j]) / length(j)
  }))
}

#' Neighborhood treatment under majority exposure
#'
#' Internal helper that computes the majority-treated indicator for each
#' observation.
#'
#' @param Di Numeric or integer vector of own treatment assignments.
#' @param neighbors List of neighbor indices for each observation.
#'
#' @return Numeric vector of majority-treatment indicators.
#'
#' @keywords internal
#' @noRd
majority <- function(Di, neighbors) {
  as.numeric(sapply(neighbors, function(j) {
    m <- length(j)
    if (m == 0L) return(NA_real_)
    sum(Di[j]) >= ceiling(m / 2)
  }))
}

#' Boundary codimension under sum-treated exposure
#'
#' Internal helper that computes the codimension of the boundary between two
#' effective treatments under the `"sum.treated"` exposure mapping.
#'
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#'
#' @return Integer codimension.
#'
#' @keywords internal
#' @noRd
ss_sum_treated <- function(d, g, dprime, gprime) {
  as.integer(d != dprime) + abs(g - gprime)
}

#' Boundary codimension under fraction-treated exposure
#'
#' Internal helper that computes the codimension of the boundary between two
#' effective treatments under the `"frac.treated"` exposure mapping.
#'
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param m Integer neighborhood size.
#' @param tol Tolerance used for checking integer compatibility.
#'
#' @return Integer codimension, or `NA_real_` if the exposure values are not
#'   compatible with neighborhood size `m`.
#'
#' @keywords internal
#' @noRd
ss_frac_treated <- function(d, g, dprime, gprime, m, tol = 1e-8) {
  k <- g * m
  kp <- gprime * m

  if (abs(k - round(k)) > tol || abs(kp - round(kp)) > tol) {
    return(NA_real_)
  }

  k <- as.integer(round(k))
  kp <- as.integer(round(kp))

  as.integer(d != dprime) + abs(k - kp)
}

#' Boundary codimension under majority exposure
#'
#' Internal helper that computes the codimension of the boundary between two
#' effective treatments under the `"majority"` exposure mapping.
#'
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param m Optional neighborhood size.
#'
#' @return Integer codimension.
#'
#' @keywords internal
#' @noRd
ss_majority <- function(d, g, dprime, gprime, m = NULL) {
  as.integer(d != dprime) + as.integer(g != gprime)
}

#' Human-readable boundary description under one-treated exposure
#'
#' Internal helper used by `get_boundary()` to generate the textual description
#' attached to a boundary object.
#'
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#'
#' @return Character scalar describing the boundary.
#'
#' @keywords internal
#' @noRd
descr_one_treated <- function(d, g, dprime, gprime) {
  dd <- d + dprime
  gg <- g + gprime

  descr <- character(2)

  if (dd == 1) descr[1] <- "Xi = 0"
  if (dd == 2) descr[1] <- "Xi >= 0"
  if (dd == 0) descr[1] <- "Xi < 0"

  if (gg == 1) descr[2] <- "there exists a j in the interference set such that Xj = 0"
  if (gg == 0) descr[2] <- "Xj < 0 for all j in the interference set."
  if (gg == 2) descr[2] <- "there exists a j in the interference set such that Xj >= 0"

  paste(descr, collapse = ", ")
}

#' Human-readable boundary description under sum-treated exposure
#'
#' Internal helper used by `get_boundary()` to generate the textual description
#' attached to a boundary object.
#'
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param dimension Integer neighborhood size.
#'
#' @return Character scalar describing the boundary.
#'
#' @keywords internal
#' @noRd
descr_sum_treated <- function(d, g, dprime, gprime, dimension) {
  dd <- d + dprime
  gg <- g - gprime
  ct <- dimension - g

  descr <- character(2)

  if (dd == 1) descr[1] <- "Xi = 0"
  if (dd == 2) descr[1] <- "Xi >= 0"
  if (dd == 0) descr[1] <- "Xi < 0"

  descr[2] <- paste0(
    "Xj = 0 for ", gg, " neighbors, ",
    "Xj >= 0 for ", gprime, " neighbors, ",
    "Xj < 0 for ", ct, " neighbors."
  )

  paste(descr, collapse = ", ")
}

#' Human-readable boundary description under fraction-treated exposure
#'
#' Internal helper used by `get_boundary()` to generate the textual description
#' attached to a boundary object.
#'
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param dimension Integer neighborhood size.
#'
#' @return Character scalar describing the boundary.
#'
#' @keywords internal
#' @noRd
descr_frac_treated <- function(d, g, dprime, gprime, dimension) {
  k <- as.integer(round(g * dimension))
  kprime <- as.integer(round(gprime * dimension))

  dd <- d + dprime
  descr <- character(2)

  if (dd == 1) descr[1] <- "Xi = 0"
  if (dd == 2) descr[1] <- "Xi >= 0"
  if (dd == 0) descr[1] <- "Xi < 0"

  if (k == kprime) {
    descr[2] <- paste0(
      "same treated-neighbor fraction (",
      k, " out of ", dimension, " treated neighbors)"
    )
  } else {
    gl <- min(k, kprime)
    gh <- max(k, kprime)
    r <- gh - gl

    descr[2] <- paste0(
      r, " neighbors at Xj = 0, ",
      gl, " neighbors with Xj >= 0, ",
      dimension - gh, " neighbors with Xj < 0"
    )
  }

  paste(descr, collapse = ", ")
}

#' Human-readable boundary description under majority exposure
#'
#' Internal helper used by `get_boundary()` to generate the textual description
#' attached to a boundary object.
#'
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param dimension Integer neighborhood size.
#'
#' @return Character scalar describing the boundary.
#'
#' @keywords internal
#' @noRd
descr_majority <- function(d, g, dprime, gprime, dimension) {
  dd <- d + dprime
  q <- ceiling(dimension / 2)

  descr <- character(2)

  if (dd == 1) descr[1] <- "Xi = 0"
  if (dd == 2) descr[1] <- "Xi >= 0"
  if (dd == 0) descr[1] <- "Xi < 0"

  if (g == gprime) {
    descr[2] <- "no change in majority exposure"
  } else if ((g == 1 && gprime == 0) || (g == 0 && gprime == 1)) {
    descr[2] <- paste0(
      "the ", q, "-th largest neighbor score equals 0 ",
      "(equivalently, exactly the threshold for a majority-treated neighborhood)"
    )
  } else {
    descr[2] <- "majority exposure must be binary with values 0/1"
  }

  paste(descr, collapse = ", ")
}
