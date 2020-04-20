#' Diagnosis plots for Lambda, U, C, Mu, Pi, Z and Epsilon
#' 
#' Diagnostics plots for hazard rate (Lambda), latent variable (U), dependence
#' variable (C), mean of cure threshold (Mu), cure proportion (Pi), cure threshold (Z) and the parameter of the
#' hierarchical prior (Epsilon).
#' 
#' This function returns a diagnostics plot for which the chain for the selected
#' variable can be monitored. Diagnostics includes trace, ergodic mean,
#' autocorrelation function and histogram.
#' 
#' @param M List. Contains the output by \code{CuMRes}.
#' @param variable Either "Lambda", "U", "C", "Mu", "Pi", "Z" or "Epsilon".
#' Variable for which diagnostic plot will be shown.
#' @param pos Positive integer. Position of the selected \code{variable} to be
#' plotted.
#' @seealso \link{CuMRes}
#' @references Nieto-Barajas, L. E., & Yin, G. (2008). Bayesian semiparametric
#' cure rate model with an unknown threshold. \emph{Scandinavian Journal of
#' Statistics}, \strong{35(3)}, 540-556.
#' https://doi.org/10.1111/j.1467-9469.2007.00589.x
#' @examples
#' 
#' 
#' 
#' ## Simulations may be time intensive. Be patient.
#' 
#' ## Example 1
#' # data(crm3)
#' # times<-crm3$times
#' # delta<-crm3$delta
#' # res <- CuMRes(times, delta, type.t = 2, 
#' #                   K = 100, length = .1, alpha = rep(1, 100  ), 
#' #                   beta = rep(1, 100),c.r = rep(50, 99), 
#' #                   iterations = 100, burn.in = 10, thinning = 1, type.c = 2)
#' # CuPlotDiag(M = res, variable = "Mu")
#' # CuPlotDiag(M = res, variable = "Z")
#' # CuPlotDiag(M = res, variable = "Pi")
#' # CuPlotDiag(M = res, variable = "Lambda", pos = 2)
#' # CuPlotDiag(M = res, variable = "U", pos = 4)
#' # CuPlotDiag(M = res, variable = "C", pos = 3)
#' 
#' 
#' 
#' @export CuPlotDiag
CuPlotDiag <-
  function(M, variable = "Lambda", pos = 1) {
    variable <-  match.arg(variable,c("Lambda","U","C","Mu","Pi","Z","Epsilon"))
    K <- extract(M, "K")
    if (pos < 0 || pos > K ) {
      stop ("Invalid position.")
    }
    if (pos > (K - 1) && (variable == "U" || variable == "C")) {
      stop ("Invalid position.")
    }
    if (!("Epsilon" %in% (dplyr::pull(extract(M, c("simulations")), name))) && variable == "Epsilon"){
      stop("Plots for 'epsilon' are not available.")
    }
    if (variable == "Epsilon" && pos != 1) {
      warning("'epsilon' has only one entry (1). Graphics shown for epsilon_1.")
      pos <- 1
    }
    MAT <- extract(M, c("simulations",variable))
    if(variable %in% c("Epsilon","Mu","Pi","Z")){
      pos = 1
      MAT <- tibble::as_tibble(matrix(MAT, nrow = extract(M, "iterations"), ncol = 1))
      
    }
    if(variable == "Lambda"){
      MAT <- rlang::set_names(dplyr::select(MAT, pos), "V1")
    } 
    if(variable %in% c("U","C")){
      MAT <- rlang::set_names(dplyr::select(tibble::as_tibble(MAT), pos), "V1")
    }
    
    var <- switch(variable, Lambda = expression(lambda),
                  Pi = expression(pi),
                  Mu = expression(mu),
                  Epsilon = expression(epsilon),
                  Z = "Z",
                  U = "U",
                  C = "C")
    if(variable %in% c("Epsilon","Mu","Pi","Z")) title <- "" else{title <- paste0("Position: ", pos)}
    a <- ggplot2::ggplot(MAT) + ggplot2::geom_line(ggplot2::aes(x=seq_len(nrow(MAT)), y = V1), color = "slateblue4") + 
      ggplot2::labs(x = "Iteration", y = variable) + ggplot2::ylab(var) + ggplot2::ggtitle("Trace")+
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    
    b <- ggplot2::ggplot(MAT) + ggplot2::geom_line(ggplot2::aes(x=seq_len(nrow(MAT)), y = cumsum(V1)/seq_len(nrow(MAT))), color = "slateblue4") + 
      ggplot2::labs(x = "Iteration", y = variable) + ggplot2::ggtitle("Ergodic mean") +
      ggplot2::ylab(var) +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    acf.aux <- acf(MAT, plot  = F)
    c <- ggplot2::ggplot(tibble::as_tibble(cbind(acf.aux$lag, acf.aux$acf))) + 
      ggplot2::geom_segment(ggplot2::aes(x = V1, xend = V1, y = V2, yend = 0)) + ggplot2::labs(x = "Lag", y ="ACF")+
      ggplot2::ggtitle("Autocorrelation function") +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    d <- ggplot2::ggplot(MAT) + ggplot2::geom_histogram(ggplot2::aes(x = V1), fill = "lightblue", color = "black", bins = 30) + 
      ggplot2::ggtitle("Histogram") + ggplot2::xlab(var) + ggplot2::ylab("") + #coord_flip() +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    gridExtra::grid.arrange(a,b,c,d, top = title)

  }
