BeUpdC <-
  function(alpha, beta, Pi.r, u.r, epsilon) {
    ck <- 50
    K <- length(Pi.r)
    c.r <- purrr::map_int(seq_len(K-1),function(index = .x){
      id <- seq.int(u.r[index], u.r[index] + ck)
      probs <- (lgamma(alpha[index + 1] + beta[index + 1] + id) 
        - lgamma(beta[index + 1] + id - u.r[index]) - lgamma(id - u.r[index] + 1)) +
        id * (log(epsilon) + log(1 - Pi.r[index]) 
              + log(1 - Pi.r[index + 1]))
      probs <- exp(probs)
      sample(x = id, size = 1, prob = probs)
    }) 
    return(c.r)
  }





