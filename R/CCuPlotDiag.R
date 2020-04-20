#' Diagnostics plots for Lambda, Theta, Delta,
#' U, C, Pi, Z and Epsilon. Hazard function, cure proportion and cure time for the median observation.
#' 
#' Diagnostic plots for hazard rate (Lambda), regression parameters for the
#' hazard (Theta), regression parameters for the cure rate (Delta), latent
#' variable (U), dependence parameter (C), mean of cure threshold (Mu), 
#' cure proportion (Pi), cure threshold (Z) and the
#' parameter of the hierarchical prior (Epsilon).
#' 
#' This function returns a diagnosyics plot for which the chain for the selected
#' variable can be monitored. Diagnostics includes trace, ergodic mean,
#' autocorrelation function and histogram.
#' 
#' @param M tibble. Contains the output by
#' \code{CCuMRes}.
#' @param variable Either "Lambda", "U", "C", "Mu", "Pi", "Z" or "Epsilon".
#' Variable for which diagnostic plot will be shown.
#' @param pos Positive integer. Position of the selected \code{variable} to be
#' plotted.
#' @seealso \link{CCuMRes}
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
#' # data(BMTKleinbook)
#' # res <- CCuMRes(BMTKleinbook, covs.x = c("tTransplant","hodgkin","karnofsky","waiting"),
#' #                covs.y = c("tTransplant","hodgkin","karnofsky","waiting"),
#' #                        type.t = 2, K = 72, length = 30,
#' #                        alpha = rep(2,72), beta = rep(2,72), c.r = rep(50, 71), type.c = 2,
#' #                        var.delta.str = .1, var.theta.str = 1,
#' #                        var.delta.ini = 100, var.theta.ini = 100,
#' #                        iterations = 100, burn.in = 10, thinning = 1)
#' # CCuPlotDiag(M = res, variable = "Z")
#' # CCuPlotDiag(M = res, variable = "Pi.m")
#' # CCuPlotDiag(M = res, variable = "Lambda", pos = 2)
#' # CCuPlotDiag(M = res, variable = "U", pos = 4)
#' 
#' 
#' 
#' 
#' @export CCuPlotDiag
CCuPlotDiag <-
  function(M, variable = "Lambda", pos = 1) {
    variable <-  match.arg(variable,c("Lambda","Lambda.m","U","C","Theta","Delta","Pi.m","Pi","Z","Z.m","Epsilon"))
    K <- extract(M, "K")
    if (pos < 0 || pos > K ) {
      stop ("Invalid position.")
    }
    if (pos > (K - 1) && (variable == "U" || variable == "C")) {
      stop ("Invalid position.")
    }
    if (pos > (K) && (variable == "Z" || variable == "Pi")) {
      stop ("Invalid observation")
    }
    if (!("Epsilon" %in% (dplyr::pull(extract(M, c("simulations")), name))) && variable == "Epsilon"){
      stop("Plots for 'epsilon' are not available.")
    }
    if (variable == "Epsilon" && pos != 1) {
      warning("'epsilon' has only one entry (1). Graphics shown for epsilon_1.")
      pos <- 1
    }
    if (variable == "Z.m" && pos != 1) {
      warning("'Z.m' has only one entry (1). Graphics shown for Z.m_1.")
      pos <- 1
    }
    if (variable == "Pi.m" && pos != 1) {
      warning("'Pi.m' has only one entry (1). Graphics shown for Pi.m_1.")
      pos <- 1
    }
    MAT <- extract(M,c("simulations",variable))
    if(variable %in% c("Lambda.m")){
      MAT <- rlang::set_names(dplyr::select(MAT[[1]],pos), "V1")
    } else{
      MAT <- rlang::set_names(dplyr::select(tibble::as_tibble(MAT),pos), "V1") 
    }
    
    
    var <- switch(variable, Lambda = expression(lambda),
                  Lambda.m = expression(lambda[median]),
                  Pi.m = expression(pi),
                  Epsilon = expression(epsilon),
                  Theta = expression(theta),
                  Delta = expression(delta),
                  Z.m = expression(Z[median]),
                  Pi = "Pi",
                  Z = "Z",
                  U = "U",
                  C = "C")
    
    title <- paste0("Position: ", pos)
    
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
      ggplot2::geom_segment(ggplot2::aes(x = V1, xend = V1, y = V2, yend = 0)) + 
      ggplot2::labs(x = "Lag", y ="ACF")+
      ggplot2::ggtitle("Autocorrelation function") +
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    d <- ggplot2::ggplot(MAT) + ggplot2::geom_histogram(ggplot2::aes(x = V1), fill = "lightblue", color = "black", bins = 30) + 
      ggplot2::ggtitle("Histogram") + ggplot2::xlab(var) + ggplot2::ylab("") + 
      ggthemes::theme_tufte() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))
    
    gridExtra::grid.arrange(a,b,c,d, top = title)
    
  }
