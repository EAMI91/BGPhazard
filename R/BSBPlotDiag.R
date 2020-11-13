#' Plot diagnostics for BSBHaz model
#'
#' @param bsbhaz An object of class 'BSBHaz' created by
#'   \code{\link{BSBHaz}}.
#' @param variable A character indicating which variable to get the plot from.
#' @param type A character indicating if the plot should be a traceplot or plot
#'   the ergodic means.
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
#' BSBPlotDiag(samples, variable = "omega1", type = "traceplot")
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr n
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove_all
#' @importFrom rlang .data
BSBPlotDiag <- function(bsbhaz,
                        variable = c("omega1", "omega2", "lambda1",
                                     "lambda2", "gamma", "theta"),
                        type = c("traceplot", "ergodic_means")) {
  
  stopifnot(inherits(bsbhaz, "BSBHaz"))
  variable <- match.arg(variable)
  type <- match.arg(type)
  p <- switch (type,
               traceplot = plot_traceplots(bsbhaz[[variable]]),
               ergodic_means = plot_ergodic_means(bsbhaz[[variable]])
  )
  
  return(p)
  
}
