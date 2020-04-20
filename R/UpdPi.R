UpdPi <-
  function(alpha, beta, c.r, u.r, n, m) {
    K <- length(alpha)
    a <- c(alpha[1] + u.r[1] + n[1],
           alpha[seq_len(K-2) + 1] + u.r[seq_len(K-2)] + u.r[seq_len(K-2) + 1] + n[seq_len(K-2) + 1],
           alpha[K] + u.r[K - 1] + n[K])
    
    b <- c(beta[1] + c.r[1] - u.r[1] + m[1],
           beta[seq_len(K-2) + 1] + c.r[seq_len(K-2)] - u.r[seq_len(K-2)] + c.r[seq_len(K-2) + 1] - u.r[seq_len(K-2) + 1] + m[seq_len(K-2) + 1],
           beta[K] + c.r[K - 1] - u.r[K - 1] + m[K])
    Pi.r <- purrr::map_dbl(rbeta(K,shape1 = a, shape2 = b),~min(max(0.001, .x), 0.999))
    return(Pi.r)
  }


