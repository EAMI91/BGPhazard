#' Diagnostics plots for lambda, U, C, Epsilon and Theta
#' 
#' Diagnostics plots for hazard rate (Lambda), latent variable (U), dependence
#' variable (C), parameter of the hierarchical model (Epsilon) and regression
#' coefficients (Theta).
#' 
#' This function returns a diagnostics plot for the chain of the selected
#' variable. The diagnostics includes trace, ergodic mean, autocorrelation
#' function and histogram.
#' 
#' @param M Tibble. Contains the output by \code{CGaMRes}
#' @param variable Either "Lambda", "U", "C", "Epsilon" or "Theta". Variable
#' for which diagnostics plot will be shown.
#' @param pos Positive integer. Position of the selected \code{variable} to be
#' plotted.
#' @seealso \link{CGaMRes}
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
#' #  data(leukemiaFZ)
#' #  leukemia1 <- leukemiaFZ
#' #  leukemia1$wbc <- log(leukemiaFZ$wbc)
#' #  CGEX1 <- CGaMRes(data = leukemia1, K = 10, iterations = 1000, thinning = 1)
#' #  CGaPlotDiag(CGEX1,variable="Theta",pos=1)
#' 
#' 
#' 
#' @export CGaPlotDiag
CGaPlotDiag <-
  function(M, variable = "Lambda", pos = 1) {
    variable <-  match.arg(variable,c("Lambda","Theta","U","C","Epsilon"))
    K <- M %>% extract("K")
    if (pos < 0 || pos > K ) {
      stop ("Invalid position.")
    }
    if (pos > (K - 1) && (variable == "U" || variable == "C")) {
      stop ("Invalid position.")
    }
    if (!("Epsilon" %in% (M %>% extract(c("simulations")) %>% 
                          dplyr::pull(name))) && variable == "Epsilon"){
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
      if(variable %in% c("Lambda")){
        MAT %<>% purrr::pluck(1) %>% tibble::as_tibble %>% dplyr::select(pos) %>% rlang::set_names("V1")
      } else{
        MAT %<>% tibble::as_tibble() %>% dplyr::select(pos) %>% rlang::set_names("V1")  
      }
      
    }
    
    var <- switch(variable, Lambda = expression(lambda),
                  Epsilon = expression(epsilon),
                  Theta = expression(theta),
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
