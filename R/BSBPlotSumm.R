#' Plot summaries for BSBHaz model
#'
#' @param bsbhaz An object of class 'BSBHaz' created by
#'   \code{\link{BSBHaz}}.
#' @param variable A character indicating the variable to plot.
#'
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
#' BSBPlotSumm(samples, "s1")
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom ggplot2 ggplot aes geom_segment theme_bw theme scale_x_continuous
#'   labs element_blank
BSBPlotSumm <- function(bsbhaz,
                        variable = c("lambda1", "lambda2", "s1", "s2")
                        ) {
  stopifnot(inherits(bsbhaz, "BSBHaz"))
  variable <- match.arg(variable)
  p <- switch (
    variable,
    lambda1 = plot_hazards(
      BSBSumm(bsbhaz, "lambda1"),
      attr(bsbhaz, "int_len"),
      "lambda1"
    ),
    lambda2 = plot_hazards(
      BSBSumm(bsbhaz, "lambda2"),
      attr(bsbhaz, "int_len"),
      "lambda2"
    ),
    s1 = plot_survival(
      BSBSumm(bsbhaz, "s1"),
      attr(bsbhaz, "int_len"),
      "S1"
    ),
    s2 = plot_survival(
      BSBSumm(bsbhaz, "s2"),
      attr(bsbhaz, "int_len"),
      "S2"
    )
  )
  
  return(p)
  
}
