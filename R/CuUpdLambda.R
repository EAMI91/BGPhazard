CuUpdLambda <-
  function(alpha, beta, c.r, u.r, n, m, z) {
    K <- length(alpha)
    tol <- 1e-7
    lambda.r <- rgamma(K,shape = c(alpha[1] + u.r[1] + n[1],
                                 alpha[seq_len(K-2) + 1] + u.r[seq_len(K-2)] + u.r[seq_len(K-2) + 1] + n[seq_len(K-2) + 1],
                                 alpha[K] + u.r[K - 1] + n[K]),
                          scale = c(1 / (beta[1] + c.r[1] + m[1]),
                                 1/(beta[seq_len(K-2) + 1] + c.r[seq_len(K-2)] + c.r[seq_len(K-2) + 1] + m[seq_len(K-2) + 1] * ((seq_len(K-2) + 1) <= z)),
                                 1 / (beta[K] + c.r[K - 1] + m[K]))
                       )
    
    lambda.r[abs(lambda.r - 0) < tol] <- 0.000001
    return(lambda.r)
  }


