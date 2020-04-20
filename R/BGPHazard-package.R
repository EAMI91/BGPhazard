#' BGPHazard: A package bayesian nonparametric inference in survival analysis.
#'
#' The BGPHazard package provides three categories of important functions:
#' simulating, diagnostic and result.
#' 
#' @section Simulating functions:
#' The simulating functions are used to make posterior inference for the bayesian survival 
#' semiparametric models as described by  
#' Nieto-Barajas and Walker (2002), Nieto-Barajas (2003) and Nieto-Barajas, L. E., & Yin, G. (2008)
#' 
#' @section Diagnostic functions:
#' The diagnostic functions are used to make convergence diagnosics plots about the simulations of the parameters/variables.
#' 
#' @section Result functions:
#' The result functions are used to produce estimators plots of the hazard function 
#' along with the survival function defined by the model.
#'
#' @docType package
#' @name BGPHazard
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang !! !!!
#' @importFrom stats acf as.formula dgamma dnorm median quantile rbeta rgamma rnorm runif time 
#' @importFrom utils data
#'

utils::globalVariables(c("name", "V1", "V2","lower","upper","S^(t)","h.est","surv",".x","value","x","y","times","."))

