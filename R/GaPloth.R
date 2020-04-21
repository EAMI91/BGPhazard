#' Plots for the Hazard and Survival Function Estimates
#' 
#' Plots the hazard function and with the survival function
#' estimates defined by the Markov gamma process with and without covariates
#' (Nieto-Barajas & Walker, 2002).
#' 
#' This function returns estimators plots for the resulting hazard rate as it is computed
#' by \link{GaMRes} and \link{CGaMRes} and the Nelson-Aalen
#' estimate along with their confidence intervals for the data set given.
#' Additionally, it plots the survival function and the Kaplan-Meier estimate
#' with their corresponding credible/confidence intervals.
#' 
#' @param M tibble. Contains the output by \code{CGaMRres} and \code{GaMRes}.
#' @param type.h character. "segment"= use segments to plot hazard rates,
#' "line" = link hazard rates by a line
#' @param addSurvival Logical. If \code{TRUE}, Nelson-Aalen estimate is plotted
#' over the hazard function and Kaplan-Meier estimate is plotted over the
#' survival function.
#' @param intervals logical. If TRUE, plots confidence bands for the selected functions including Nelson-Aalen and/or Kaplan-Meier estimate.
#' @param confidence Numeric. Confidence level.
#' @param summary Logical. If \code{TRUE}, a summary for hazard and survival
#' functions is returned as a tibble.
#' @return \item{SUM.h}{Numeric tibble. Summary for the mean, median, and a
#' \code{confint / 100} confidence interval for each segment of the hazard
#' function. If \code{summary = TRUE}} \item{SUM.S}{Numeric tibble. Summary for
#' the mean, median, and a \code{confint / 100} confidence interval for a grid
#' of the survival function. If \code{summary = TRUE}}
#' @seealso \link{GaMRes}, \link{CGaMRes}, \link{CGaPlotDiag},
#' \link{GaPlotDiag}
#' @references - Nieto-Barajas, L. E. (2003). Discrete time Markov gamma
#' processes and time dependent covariates in survival analysis. \emph{Bulletin
#' of the International Statistical Institute 54th Session}. Berlin. (CD-ROM).
#' 
#' - Nieto-Barajas, L. E. & Walker, S. G. (2002). Markov beta and gamma
#' processes for modelling hazard rates. \emph{Scandinavian Journal of
#' Statistics} \strong{29}: 413-424.
#' @examples
#' 
#' 
#' 
#' ## Simulations may be time intensive. Be patient.
#' 
#' ## Example 1
#' #  data(gehan)
#' #  timesG <- gehan$time[gehan$treat == "6-MP"]
#' #  deltaG <- gehan$cens[gehan$treat == "6-MP"]
#' #  GEX1 <- GaMRes(timesG, deltaG, K = 8, iterations = 3000)
#' #  GaPloth(GEX1)
#' 
#' 
#' ## Example 2
#' #  data(leukemiaFZ)
#' #  timesFZ <- leukemiaFZ$time
#' #  deltaFZ <- leukemiaFZ$delta
#' #  GEX2 <- GaMRes(timesFZ, deltaFZ, type.c = 4)
#' #  GaPloth(GEX2)
#' 
#' 
#' 
#' 
#' 
#' @export GaPloth
GaPloth <-
  function(M, type.h = "segment", addSurvival = T, intervals = T,
           confidence = 0.95, summary = FALSE) {
    SUM <- LambdaSumm(M, confidence)
    s <- SUM %>% tibble::deframe()
    v <- list("tao",
              "K",
              "times",
              "delta"
    ) %>% purrr::map(~extract(M,.x)) %>% rlang::set_names(c("tao","K","times","delta"))
    tao <- v$tao
    K <- v$K
    delta <- v$delta
    times <- v$times
    
    if(type.h == "segment") {
      h <- s$SUM.h %>% ggplot2::ggplot() + 
      ggplot2::geom_segment(ggplot2::aes(x = tao[-(K+1)], xend = tao[-1], 
                       y = mean, yend = mean, color = "Hazard Function")) + 
      ggplot2::scale_color_manual(values = c("black"), limits = "Hazard Function") +
      ggplot2::guides(color = ggplot2::guide_legend(title = "")) +
      ggplot2::xlab("Time") +ggplot2::ylab("Hazard rate") + ggplot2::scale_alpha_continuous(guide = F) + 
      ggplot2::ggtitle(paste0("Estimate of hazard rates with intervals at ",confidence * 100,"% of credibility")) +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
            legend.position="bottom")
      if(intervals){
        h <- h + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper, x = (tao[-(K+1)] + tao[-1])/2, width = tao[-1]-tao[-(K+1)]), 
                                        alpha = 0.5, color = "gray50")
      }
    }
    
    if(type.h == "line"){ 
      h <- s$SUM.h %>% ggplot2::ggplot() + 
      ggplot2::geom_line(ggplot2::aes(x = (tao[-(K+1)] + tao[-1])/2, y = mean, color = "Hazard Function")) +
      ggplot2::scale_color_manual(values = c("black"), limits = "Hazard Function") +
      ggplot2::guides(color = ggplot2::guide_legend(title = "")) +
      ggplot2::xlab("Time") + ggplot2::ylab("Hazard rate") + ggplot2::scale_alpha_continuous(guide = F) + 
      ggplot2::ggtitle(paste0("Estimate of hazard rates with intervals at ",confidence * 100,"%  of credibility")) +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
            legend.position="bottom")
      if(intervals){
        h <- h + ggplot2::geom_ribbon(ggplot2::aes(x = (tao[-(K+1)] + tao[-1])/2, ymin = lower, ymax = upper), alpha = .5, fill = "gray70")
      }
    }
    
    S <- s$SUM.S %>% ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x = t, y = `S^(t)`,color = "Model estimate")) + 
      ggplot2::scale_color_manual(limits = c("Model estimate"),values = c("black")) +
      ggplot2::guides(color = ggplot2::guide_legend(title = "")) +
      ggplot2::scale_y_continuous(limits = c(0,1)) + 
      ggplot2::ggtitle(paste0("Estimate of Survival Function with intervals at ", confidence * 100,"%  of credibility")) +
      ggplot2::labs(x = "t",
           y = expression(S^{(t)})) +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
            legend.position = "bottom")
    
    if(intervals){
      S <- S + ggplot2::geom_ribbon(ggplot2::aes(x = t, ymin = lower, ymax = upper), fill = "gray50", alpha = 0.3)
    }
    
    if(addSurvival){
      fit <- survival::survfit(survival::Surv(time = times, event = delta) ~ 1,
                     conf.int = confidence)
      
      km.data <- tibble::tibble(time = fit$time,surv = fit$surv, lower = fit$lower,
                        upper = fit$upper)
      if(km.data$time[1]!= 0){
        km.data <- dplyr::bind_rows(tibble::tibble(time = 0, surv = 1, lower = NA, upper = NA),km.data)
      }
      na.data <- tibble::tibble(time = fit$time, h.est = fit$n.event / fit$n.risk)
      h <- h + ggplot2::geom_point(data = na.data, ggplot2::aes(x = time, y = h.est), color = "#b22222") + 
        ggplot2::scale_color_manual(limits = c("Hazard Function","Nelson-Aalen based estimate"),
                           values = c("black","#b22222")) 
      
      S <- S + ggplot2::geom_step(data = km.data,na.rm = T, ggplot2::aes(x = time,y = surv), color = "#b22222") + 
        ggplot2::scale_color_manual(limits = c("Model estimate","Kaplan Meier"),
                           values = c("black","#b22222")) 
      
      if(intervals){
        S <- S + ggplot2::geom_step(data = km.data, na.rm = T, ggplot2::aes(x = time, y = lower), alpha = 0.5, color = "#b22222", linetype = "dashed") + 
          ggplot2::geom_step(data = km.data, na.rm = T, ggplot2::aes(x = time, y = upper), alpha = 0.5, color = "#b22222", linetype = "dashed")
      }
    }
    
    if (summary == TRUE) {
      return(list(h,S,SUM))
    } else{
      return(list(h,S))
    }
  }
