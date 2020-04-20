BeUpdU <-
  function(alpha,  beta,  c.r,  Pi.r) {
    
    ind <- which(c.r != 0)
    u.r <- rep(0,length(c.r))
    if(length(ind)>0){
      K <- length(Pi.r)
      lphi <- (log(Pi.r[-K]) + log(Pi.r[-1]) 
               - log(1 - Pi.r[-K]) - log(1 - Pi.r[-1]))
      
      aux_u <- purrr::map(ind,
          ~exp(seq.int(0, c.r[.x]) * lphi[.x] - 
                 lgamma(seq.int(0, c.r[.x]) + 1) - 
                 lgamma(c.r[.x] - seq.int(0, c.r[.x]) + 1) -
                 lgamma(alpha[.x + 1] + seq.int(0,c.r[.x])) -
                 lgamma(beta[.x + 1] + c.r[.x] - seq.int(0,c.r[.x])))
      )
      u.r[ind] <-  purrr::map_int(aux_u, 
        ~sample(x = seq.int(0,length(.x)-1),size = 1, prob = .x))  
    }
    
    return(u.r)
  }






