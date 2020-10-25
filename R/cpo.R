#' Conditional Predictive Ordinate (CPO) Statistic
#' 
#' Makes the CPO Plot and calculates the logarithm of the Pseudomarginal
#' likelihood (LPML).
#' 
#' Computes de CPO as a goodness of fit measure
#' 
#' @param res tibble. The output from the *Res functions, where * could either
#' be BeM, GaM, CGaM, CuM, CCuM
#' @return 
#' \item{LPML}{The value of the logarithm of the Pseudomarginal likelihood}
#' \item{plot}{CPO Plot} %% ...
#' @references See Geisser (1993); Gelfand, Dey, and Chang (1992); Dey, Chen,
#' and Chang (1997); and Sinha and Dey (1997)
#' @examples
#' 
#' 
#' 
#' ## Example 1
#' #  data(gehan)
#' #  timesG <- gehan$time[gehan$treat == "6-MP"]
#' #  deltaG <- gehan$cens[gehan$treat == "6-MP"]
#' #  GEX1 <- GaMRes(timesG, deltaG, K = 8, iterations = 3000)
#' #  cpo(GEX1)
#' 
#' 
#' 
#' @export cpo
cpo <- function(res){
  #Extraer variables necesarias para el cálculo de CPO
  aux <- res %>% extract("simulations") %>% dplyr::pull(name)
  if("Lambda" %in% aux) parameter <- "Lambda" else parameter <- "PI"
  aux <- list("times",
              "delta",
              "K",
              "tao",
              c("simulations",parameter),
              "s",
              "S"
  ) %>% purrr::map(~res %>% extract(.x)) %>% rlang::set_names("times","delta","k","tao","lambda","s","S")
  
  if(length(aux$S) == 1) aux$S <- aux$S %>% purrr::pluck(1)
  if(length(aux$lambda) == 1) aux$lambda <- aux$lambda %>% purrr::pluck(1)
  uncensored <- aux$times[aux$delta == 1] %>%
    cut(aux$tao,labels = F,include.lowest = T) %>%
    purrr::map2_dbl(.y = aux$times[aux$delta == 1] %>%
               cut(aux$s,labels = F,include.lowest = T), 
             ~(mean(((tibble::as_tibble(aux$lambda) %>% dplyr::pull(.x)) * (aux$S %>% purrr::pluck(.y)))^(-1)))^(-1)
    )
  
  censored <- aux$times[aux$delta == 0] %>%
    cut(aux$s,labels = F,include.lowest = T) %>% 
    purrr::map_dbl(~mean(aux$S[[.x]]^(-1))^(-1) ) 
  
  g <- tibble::tibble(times=c(aux$times[aux$delta == 1],
                      aux$times[aux$delta == 0]), cpo = c(uncensored,censored)) %>%
    ggplot2::ggplot(ggplot2::aes(x=times,y=cpo))+ ggplot2::geom_point() 
  
  lpml <- c(uncensored,censored) %>% log %>% sum 
  return(list(LPML = lpml,plot = g))
}
