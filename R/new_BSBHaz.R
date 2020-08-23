# Not exported
new_BSBHaz <- function(l = list(),
                       individuals = integer(),
                       intervals = integer(),
                       has_predictors = logical(),
                       samples = integer(),
                       int_len = double()) {
  
  stopifnot(is.list(l))
  stopifnot(is.integer(individuals))
  stopifnot(is.integer(intervals))
  stopifnot(is.logical(has_predictors))
  stopifnot(is.integer(samples))
  
  structure(l,
            class = "BSBHaz",
            individuals = individuals,
            intervals = intervals,
            has_predictors = has_predictors,
            samples = samples,
            int_len = int_len
  )
  
}
