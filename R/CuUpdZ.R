CuUpdZ <-
  function(mu, m, lambda.r, k.star) {
    k <- length(lambda.r)
    propfz <- purrr::map_dbl(seq.int(k.star,k), function(z){
      logfz <- (z - 1) * log(mu) - lgamma(z) - sum( m[seq_len(z)] * lambda.r[seq_len(z)] )
      fz <- exp(logfz)
      return(fz)
    })
    if(length(propfz) == 1) z <- k else{
      z <- sample(x = seq.int(k.star,k), size = 1, prob = propfz)
    }
    return(z)
  }