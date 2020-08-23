#' @export
print.BSBinit <- function(x, ...) {
  stopifnot(inherits(x, "BSBinit"))
  individuals <- attr(x, "individuals")
  intervals <- attr(x, "intervals")
  has_predictors <- attr(x, "has_predictors")
  max_t <-
    if (has_predictors) {
      pred_names <- colnames(x$pred_matrix)
    } else {
      pred_names <- " "
    }

  cat(
    "\t\n",
    sprintf("Individuals: %s\n", individuals),
    sprintf("Time partition intervals: %s\n", intervals),
    sprintf("Censored t1: %s\n", individuals - sum(x$delta1)),
    sprintf("Censored t2: %s\n", individuals - sum(x$delta2)),
    "Predictors:", as.character(has_predictors), "\t", pred_names
  )

}

#' @export
print.BSBHaz <- function(x, ...) {
  stopifnot(inherits(x, "BSBHaz"))
  individuals <- attr(x, "individuals")
  intervals <- attr(x, "intervals")
  has_predictors <- attr(x, "has_predictors")
  samples <- attr(x, "samples")
  max_t <-
    if (has_predictors) {
      pred_names <- rownames(x$theta_mat)
    } else {
      pred_names <- " "
    }

  cat(
    "\t\n",
    sprintf("Samples: %s\n", samples),
    sprintf("Individuals: %s\n", individuals),
    sprintf("Time partition intervals: %s\n", intervals),
    "Predictors:", as.character(has_predictors), "\t", pred_names
  )
}
