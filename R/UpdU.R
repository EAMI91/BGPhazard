UpdU <-
  function(alpha, beta, c.r, lambda.r) {
    uk <- 100
    ind <- which(c.r != 0)
    u <- rep(0,length(c.r))
    if(length(ind)>0){
      aux_u <- purrr::map(ind,function(k) {
        exp(seq.int(0,uk) * (log(c.r[k]) + log(c.r[k] + beta[k + 1]) + 
                            log(lambda.r[k]) + log(lambda.r[k + 1])) - 
           (lgamma(seq.int(0,uk)+1) + lgamma(alpha[k + 1] + seq.int(0,uk)))) 
      })
      u[ind] <- purrr::map_dbl(.x = aux_u, .f=~sample(x=seq.int(0,uk), size=1, prob=.x)) 
    }
    return(u)
  }