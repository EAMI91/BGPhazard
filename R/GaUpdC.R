GaUpdC <-
  function(alpha, beta, c.r, lambda.r, u.r, epsilon, nu, acceptance.c) {
    K <- length(lambda.r)
    
    c.str <- rgamma(K-1, shape = nu, 
                    rate = nu/c.r)
    
    lw.c.str <- (alpha[seq_len(K-1) + 1] + u.r[seq_len(K-1)]) * log(beta[seq_len(K-1) + 1] + c.str) +
      u.r[seq_len(K-1)]*log(c.str) - (lambda.r[seq_len(K-1) + 1] + lambda.r[seq_len(K-1) ] + 1/epsilon)*c.str + 
      log(dgamma(x = c.r, shape = nu, rate = nu/c.str))
    lw.c.r <- (alpha[seq_len(K-1) + 1] + u.r[seq_len(K-1)]) * log(beta[seq_len(K-1) + 1] + c.r) +
      u.r[seq_len(K-1)]*log(c.r) - (lambda.r[seq_len(K-1) + 1] + lambda.r[seq_len(K-1) ] + 1/epsilon)*c.r + 
      log(dgamma(x = c.str, shape = nu, rate = nu/c.r))
    ratio <- lw.c.str - lw.c.r
    unifs <- runif(K-1)
    criteria <- log(unifs) <= ratio
    c.r[criteria] <- c.str[criteria]
    acceptance.c <- acceptance.c + sum(criteria)
    return(list(c.r, acceptance.c))
  }
