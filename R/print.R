#' Print rddinterf estimation results
#'
#' Displays the estimation method and design settings followed by a compact
#' table of effect estimates. The full boundary objects, distance variables,
#' bandwidths, and estimation samples remain available in the fitted object.
#'
#' @param x An object of class `"rddinterf"`.
#' @param digits Number of digits shown for estimates, standard errors, and
#'   confidence intervals.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.rddinterf <- function(x, digits = 3, ...) {
  if (!is.numeric(digits) || length(digits) != 1L ||
      !is.finite(digits) || digits < 0) {
    stop("digits must be one nonnegative finite number.")
  }
  digits <- as.integer(digits)

  estimator <- x$estimator
  if (is.null(estimator)) {
    estimator <- if (is.null(x$call$estimator)) {
      "conventional"
    } else {
      as.character(x$call$estimator)
    }
  }
  estimator_label <- switch(
    estimator,
    conventional = "Conventional local polynomial",
    robust = "Robust bias-corrected local polynomial",
    estimator
  )

  cat("Regression discontinuity estimates under interference\n")
  cat("Estimator: ", estimator_label, "\n", sep = "")
  cat("Exposure mapping: ", x$exposure, "\n", sep = "")
  cat("Cutoff: ", format(x$cutoff), "\n", sep = "")
  cat("Polynomial orders: p = ", x$p, ", q = ", x$q, "\n\n", sep = "")

  table_data <- summary(x)
  if (nrow(table_data) == 0L) {
    cat("No estimable effects.\n")
    return(invisible(x))
  }

  display <- data.frame(
    Estimate = formatC(table_data[, "Estimate"], format = "f", digits = digits),
    `Std. Error` = formatC(table_data[, "S.E."], format = "f", digits = digits),
    `95% CI` = paste0(
      "[",
      formatC(table_data[, "CI.l"], format = "f", digits = digits),
      ", ",
      formatC(table_data[, "CI.u"], format = "f", digits = digits),
      "]"
    ),
    `p-value` = format.pval(
      table_data[, "p.value"],
      digits = max(1L, digits),
      eps = 10^(-max(1L, digits))
    ),
    `N left` = as.integer(table_data[, "N.l"]),
    `N right` = as.integer(table_data[, "N.r"]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(display) <- rownames(table_data)

  print(display, quote = FALSE, right = TRUE)
  invisible(x)
}
