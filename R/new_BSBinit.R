# Not exported

# Creates an object of class 'BSBinit'
new_BSBinit <- function(l = list(),
                        individuals = integer(),
                        intervals = integer(),
                        has_predictors = logical()) {

  stopifnot(is.list(l))
  stopifnot(is.integer(individuals))
  stopifnot(is.integer(intervals))
  stopifnot(is.logical(has_predictors))

  structure(
    l,
    class = "BSBinit",
    individuals = individuals,
    intervals = intervals,
    has_predictors = has_predictors
    )
}
