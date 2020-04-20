#' Diagnosis plots for Lambda, U, C and Epsilon
#' 
#' Diagnostics plots for hazard rate (Lambda), latent variable (U), dependence
#' parameter (C) and the parameter of the hierarchical prior (Epsilon).
#' 
#' This function returns a diagnostics plot for which the chain of the selected
#' variable can be monitored. Diagnostics includes trace, ergodic mean,
#' autocorrelation function and histogram.
#' 
#' @param M List. Contains the output
#' by \code{GaMRes}.
#' @param variable Either "Lambda", "U", "C" or "Epsilon". Variable for which
#' informative plot will be shown.
#' @param pos Positive integer. Position of the selected \code{variable} to be
#' plotted.
#' @seealso \link{GaMRes}
#' @references - Nieto-Barajas, L. E. & Walker, S. G. (2002). Markov beta and
#' gamma processes for modelling hazard rates. \emph{Scandinavian Journal of
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
#' #  GaPlotDiag(GEX1, variable = "Lambda", pos = 2)
#' #  GaPlotDiag(GEX1, variable = "U", pos = 5)
#' 
#' ## Example 2
#' #  data(leukemiaFZ)
#' #  timesFZ <- leukemiaFZ$time
#' #  deltaFZ <- leukemiaFZ$delta
#' #  GEX2 <- GaMRes(timesFZ, deltaFZ, type.c = 4)
#' #  GaPlotDiag(GEX2, variable = "Lambda", pos = 2)
#' #  GaPlotDiag(GEX2, variable = "U", pos = 3)
#' 
#' 
#' 
#' @export GaPlotDiag
GaPlotDiag <-
  function(M, variable = "Lambda", pos = 1) {
    variable <-  match.arg(variable,c("Lambda","U","C","Epsilon"))
    K <- extract(M,"K")
    if (pos < 0 || pos > K ) {
      stop ("Invalid position.")
    }
    if (pos > (K - 1) && (variable == "U" || variable == "C")) {
      stop ("Invalid position.")
    }
    if (!("Epsilon" %in% (M %>% extract(c("simulations")) %>% dplyr::pull(name))) && variable == "Epsilon"){
      stop("Plots for 'epsilon' are not available.")
    }
    if (variable == "Epsilon" && pos != 1) {
      warning("'epsilon' has only one entry (1). Graphics shown for epsilon_1.")
      pos <- 1
    }
    MAT <- M %>% extract(c("simulations",variable))
    if(variable %in% c("Epsilon")){
      pos = 1
      MAT <- matrix(MAT, nrow = M %>% extract("iterations"), ncol = 1) %>% tibble::as_tibble()
    } else{
      MAT %<>% tibble::as_tibble() %>% dplyr::select(pos) %>% rlang::set_names("V1")
    }
    
    var <- switch(variable, Lambda = expression(lambda),
                  Epsilon = expression(epsilon),
                  U = "U",
                  C = "C")
    if(variable %in% c("Epsilon")) title <- "" else{title <- paste0("Position: ", pos)}
    a <- MAT %>% ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=seq_len(nrow(MAT)), y = V1), color = "slateblue4") + 
      ggplot2::labs(x = "Iteration", y = variable) + ggplot2::ylab(var) + ggplot2::ggtitle("Trace")+
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    
    b <- MAT %>% ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=seq_len(nrow(MAT)), y = cumsum(V1)/seq_len(nrow(MAT))), color = "slateblue4") + 
      ggplot2::labs(x = "Iteration", y = variable) + ggplot2::ggtitle("Ergodic mean") +
      ggplot2::ylab(var) +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    acf.aux <- acf(MAT, plot  = F)
    c <- cbind(acf.aux$lag, acf.aux$acf) %>% tibble::as_tibble() %>% ggplot2::ggplot() + 
      ggplot2::geom_segment(ggplot2::aes(x = V1, xend = V1, y = V2, yend = 0)) + ggplot2::labs(x = "Lag", y ="ACF")+
      ggplot2::ggtitle("Autocorrelation function") +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    d <- MAT %>% ggplot2::ggplot() + ggplot2::geom_histogram(ggplot2::aes(x = V1), fill = "lightblue", color = "black", bins = 30) + 
      ggplot2::ggtitle("Histogram") + ggplot2::xlab(var) + ggplot2::ylab("") + #coord_flip() +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    gridExtra::grid.arrange(a,b,c,d, top = title)
    

  }
