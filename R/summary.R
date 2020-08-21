#' @export
summary.BSBinit <- function(object, ...) {
  stopifnot(inherits(object, "BSBinit"))
  individuals <- attr(object, "individuals")
  intervals <- attr(object, "intervals")
  has_predictors <- attr(object, "has_predictors")
  max_t <-
  if (has_predictors) {
    pred_names <- colnames(object$pred_matrix)
  } else {
      pred_names <- " "
    }

  cat(
    "\t\n",
    sprintf("Individuals: %s\n", individuals),
    sprintf("Time partition intervals: %s\n", intervals),
    sprintf("Censored t1: %s\n", individuals - sum(object$delta1)),
    sprintf("Censored t2: %s\n", individuals - sum(object$delta2)),
    "Predictors:", as.character(has_predictors), "\t", pred_names
  )

}
