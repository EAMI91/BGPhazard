#' Get posterior summaries for BSBHaz model
#'
#' @param bsbhaz An object of class 'BSBHaz' created by
#'   \code{\link{BSBHaz}}.
#' @param variable A character indicating which variable to get summaries from.
#'
#' @return A data frame with posterior sample means and a 95 \% probability
#'   interval. For \code{omega1}, \code{omega2}, \code{gamma}, and \code{theta}
#'   also includes a column with the acceptance rates for the
#'   Metropolis-Hastings algorithm.
#' @export
#'
#' @examples
#' t1 <- survival::Surv(c(1, 2, 3))
#' t2 <- survival::Surv(c(1, 2, 3))
#'
#' init <- BSBInit(t1 = t1, t2 = t2, seed = 0)
#' samples <- BSBHaz(init, iter = 10, omega_d = 2,
#' gamma_d = 10, seed = 10)
#'
#' BSBSumm(samples, variable = "gamma")
#' BSBSumm(samples, variable = "omega1")
#' BSBSumm(samples, variable = "lambda1")
BSBSumm <- function(bsbhaz,
                    variable = c("omega1", "omega2", "lambda1",
                                 "lambda2", "gamma", "theta",
                                 "s1", "s2")) {
  stopifnot(inherits(bsbhaz, "BSBHaz"))
  variable <- match.arg(variable)
  switch (variable,
          omega1 = summaries_omega(bsbhaz$omega1),
          omega2 = summaries_omega(bsbhaz$omega2),
          lambda1 = summaries_lambda(bsbhaz$lambda1),
          lambda2 = summaries_lambda(bsbhaz$lambda2),
          theta = summaries_theta(bsbhaz$theta),
          gamma = summaries_gamma(bsbhaz$gamma),
          s1 = summaries_surv(bsbhaz$s1),
          s2 = summaries_surv(bsbhaz$s2)
  )
  
}
