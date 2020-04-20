#' Diagnosis plots for PI, U, C and Epsilon
#' 
#' Diagnostic plots for hazard rate (PI), latent variable (U), dependence
#' parameter (C) and parameter of the hierarchical model (Epsilon).
#' 
#' This function returns a diagnostics plot for the chain of the selected
#' variable. The diagnostics includes trace, ergodic mean, autocorrelation
#' function and histogram.
#' 
#' @param M Tibble. Contains the output by
#' \code{BeMRes}
#' @param variable Either "PI", "U", "C" or "Epsilon". Variable for which
#' diagnostic plot will be shown.
#' @param pos Positive integer. Position of the selected \code{variable} to be
#' plotted.
#' @seealso \link{BeMRes}
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
#' #  data(psych)
#' #  timesP <- psych$time
#' #  deltaP <- psych$death
#' #  BEX1 <- BeMRes(timesP, deltaP, iterations = 3000, burn.in = 300, thinning = 1)
#' #  BePlotDiag(BEX1, variable = "PI", pos = 2)
#' #  BePlotDiag(BEX1, variable = "U", pos = 3)
#' 
#' ## Example 2
#' #  data(gehan)
#' #  timesG <- gehan$time[gehan$treat == "control"]
#' #  deltaG <- gehan$cens[gehan$treat == "control"]
#' #  BEX2 <- BeMRes(timesG, deltaG, type.c = 2, c.r = rep(50, 22))
#' #  BePlotDiag(BEX2, variable = "PI", pos = 5)
#' #  BePlotDiag(BEX2, variable = "U", pos = 4)
#' 
#' 
#' 
#' @export BePlotDiag
BePlotDiag <-
  function(M, variable = "PI", pos = 1) {
    variable <-  match.arg(variable,c("U","C","PI","Epsilon"))
    
    K <- M %>% extract("K")
    tol = .Machine$double.eps ^ 0.5
    if (pos < 0 || pos > K || abs(pos - round(pos)) > tol ) {
      stop ("Invalid position.")
    }
    if ((variable == "C" || variable == "U") && pos > K - 1) {
      stop ("Invalid position.")
    }
    if (!("Epsilon" %in% (M %>% extract("simulations") %>% dplyr::pull(name))) && variable == "Epsilon"){
      stop("Plots for 'Epsilon' are not available.")
    }
    if (variable == "Epsilon" && pos != 1) {
      warning("'Epsilon' has only one entry (1). Graphics shown for Epsilon_1.")
      pos <- 1
    }
    MAT <- M %>% extract(c("simulations",variable))
    if(variable %in% c("Epsilon")){
      pos = 1
      MAT <- matrix(MAT, nrow = M %>% extract("iterations"), ncol = 1) %>% tibble::as_tibble()
      
    }

    if(variable %in% c("PI","U","C")){
      MAT %<>% tibble::as_tibble() %>% dplyr::select(pos) %>% rlang::set_names("V1")
    }
    
    var <- switch(variable, 
                  PI = expression(pi),
                  Epsilon = expression(epsilon),
                  U = "U",
                  C = "C")
    if(variable %in% c("Epsilon")) title <- "" else{title <- paste0("Position: ", pos)}
    a <- MAT %>% ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=seq_len(nrow(MAT)), y = V1), color = "slateblue4") + 
      ggplot2::labs(x = "Iteration", y = variable) + ggplot2::ylab(var) + ggplot2::ggtitle("Trace") +
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
