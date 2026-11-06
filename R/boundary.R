#' Construct effective treatment boundaries
#'
#' Builds boundary objects describing the separating set between two effective
#' treatment regions under a specified exposure mapping.
#'
#' @param neighbors List of neighbor indices for each observation.
#' @param contrast A length-2 list describing the two effective treatments to
#'   compare.
#' @param exposure Character string specifying the exposure mapping. Supported
#'   values currently include `"one.treated"`, `"sum.treated"`,
#'   `"frac.treated"`, and `"majority"`.
#'
#' @return An object of class `"rdd_boundary"`. Each element corresponds to a
#'   distinct neighborhood size and contains the information needed to describe
#'   the effective treatment boundary.
#'
#' @examples
#' neighbors <- list(c(2, 3), c(1, 3), c(1, 2))
#' contrast <- list(c(1, 1), c(0, 0))
#' get_boundary(neighbors, contrast, exposure = "one.treated")
#'
#' @export
get_boundary <- function(neighbors, contrast, exposure) {
  d <- contrast[[1]][1]
  g <- contrast[[1]][2]
  dprime <- contrast[[2]][1]
  gprime <- contrast[[2]][2]

  dimensions <- sort(unique(vapply(neighbors, length, integer(1))))
  boundary <- vector("list", length(dimensions))

  for (i in seq_along(dimensions)) {
    m <- dimensions[i]

    tmp <- make_boundary_entry(
      exposure = exposure,
      d = d,
      g = g,
      dprime = dprime,
      gprime = gprime,
      m = m
    )

    boundary[[i]] <- tmp$boundary

    attr(boundary[[i]], "description") <- paste(
      paste0("Effective treatment boundary of dimension: ", m + 1 - tmp$ss),
      paste0("Codimension: ", tmp$ss),
      paste0("Type of exposure: ", exposure),
      paste0("Description: ", tmp$description),
      sep = "\n"
    )
  }

  names(boundary) <- paste(dimensions, "-dimensional")
  class(boundary) <- "rdd_boundary"
  boundary
}

#' Build one boundary entry for a fixed neighborhood size
#'
#' Internal helper used by `get_boundary()` to construct the boundary metadata
#' for a given exposure mapping and neighborhood size.
#'
#' @param exposure Character string specifying the exposure mapping.
#' @param d,g Treated effective treatment components.
#' @param dprime,gprime Control effective treatment components.
#' @param m Integer neighborhood size.
#'
#' @return A list with components:
#' - `boundary`: one boundary entry in package format.
#' - `ss`: codimension of the boundary.
#' - `description`: human-readable description of the boundary.
#'
#' @keywords internal
#' @noRd
make_boundary_entry <- function(exposure, d, g, dprime, gprime, m) {
  if (exposure == "one.treated") {
    ss <- as.integer(d != dprime) + as.integer(g != gprime)
    ds <- descr_one_treated(d, g, dprime, gprime)

    entry <- list(
      exposure = exposure,
      d = d,
      g = g,
      dprime = dprime,
      gprime = gprime,
      dim = m,
      ss = ss
    )
  } else if (exposure == "sum.treated") {
    ss <- ss_sum_treated(d, g, dprime, gprime)
    ds <- descr_sum_treated(d, g, dprime, gprime, m)

    entry <- list(
      exposure = exposure,
      d = d,
      g = g,
      dprime = dprime,
      gprime = gprime,
      dim = m,
      ss = ss
    )
  } else if (exposure == "frac.treated") {
    ss <- ss_frac_treated(d, g, dprime, gprime, m)
    ds <- descr_frac_treated(d, g, dprime, gprime, m)

    entry <- list(
      exposure = exposure,
      d = d,
      g = g,
      dprime = dprime,
      gprime = gprime,
      dim = m,
      k = as.integer(round(g * m)),
      kprime = as.integer(round(gprime * m)),
      ss = ss
    )
  } else if (exposure == "majority") {
    if (!(g %in% c(0, 1)) || !(gprime %in% c(0, 1))) {
      stop("For majority exposure, g and gprime must be 0 or 1.")
    }

    ss <- ss_majority(d, g, dprime, gprime, m)
    ds <- descr_majority(d, g, dprime, gprime, m)

    entry <- list(
      exposure = exposure,
      d = d,
      g = g,
      dprime = dprime,
      gprime = gprime,
      dim = m,
      q = ceiling(m / 2),
      ss = ss
    )
  } else {
    stop("Unsupported exposure.")
  }

  list(
    boundary = list(entry),
    ss = ss,
    description = ds
  )
}


#' Print a boundary object
#'
#' Displays a human-readable representation of a boundary object, listing
#' each interference set along with its size and description.
#'
#' Boundary objects are typically created using \code{\link{get_boundary}}.
#'
#' @param x An object of class \code{"rdd_boundary"}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.rdd_boundary <- function(x, ...) {
  for (i in seq_along(x)) {
    nm <- as.numeric(regmatches(names(x)[i], regexpr("\\d+", names(x)[i])))
    if (is.na(nm)) nm <- i

    cat("\n[[ Interference set with size ", nm, "]]\n", sep = "")

    d <- attr(x[[i]], "description")
    if (!is.null(d)) {
      cat(d, "\n\n")
    }
  }
  invisible(x)
}
