CCuW <- function(theta, times, K, covar, tao, ind) {
  w <- purrr::map2(times, .y = seq_len(ind),
                 .f= ~ (.x > tao[-1]) * as.vector(exp(sum(theta * covar[.y, ]))) * (tao[-1]) + 
                   (.x > tao[-(K+1)] & .x <= tao[-1]) * as.vector(exp(theta %*% covar[.y, ])) * .x -
                   (.x > tao[-(K+1)]) * as.vector(exp(theta %*% covar[.y, ])) * tao[-(K+1)]) 
  return(w)
}
