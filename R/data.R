#' Recurrent infection of kidney catheters
#'
#' Data on the recurrent times to infection, at the point of insertion of the
#' catheter, for kidney patients using portable dialysis equipment. Catheters
#' may be removed for reasons other than infection, in which case the
#' observation is censored. Each patient has exactly 2 observations. Only sex
#' was kept as an explanatory variable.
#'
#' @format A data frame with 38 rows and 6 variables:
#' \describe{
#'   \item{id}{patient ID}
#'   \item{t1,t2}{times to infection}
#'   \item{delta1,delta2}{censorship indicators (1 = exact, 0 = right-censored)}
#'   \item{sex}{0 = female, 1 = male}
#'   }
#' @source \url{https://www.mayo.edu/research/documents/kidneyhtml/doc-10027569}
"KIDNEY"
