CCuUpdTheta <-
  function(theta, lambda.r, times, delta, K, covar, tao, ind, z, var.theta.str, var.theta.ini, acceptance.th) {
    p <- length(theta)
    theta.upd <- theta
    m.upd <- purrr::reduce(purrr::map2(CCuW(theta.upd, times, K, covar, tao, ind), .y = z, ~(seq_len(K) <= .y)*.x), `+`)
    for(s in 1:p){
      theta.str <- theta.upd
      theta.str[s] <- rnorm(1, mean = theta.upd[s], sd = sqrt(var.theta.str))
      m.str <- purrr::reduce(purrr::map2(CCuW(theta.str, times, K, covar, tao, ind), 
                                         .y = z, ~(seq_len(K) <= .y)*.x), `+`)
      pr <- dnorm(theta.str[s], mean = 0, sd = sqrt(var.theta.ini), log = T) -
        dnorm(theta.upd[s], mean = 0, sd = sqrt(var.theta.ini), log = T) +
        (theta.str[s] - theta.upd[s]) * sum(covar[delta==1, s]) + 
        sum(lambda.r * (m.upd - m.str))
      if(log(runif(1)) <= pr) {
        theta.upd <- theta.str
        m.upd <- m.str
        acceptance.th[s] <- acceptance.th[s] + 1
      }
    }
    return(list(theta.upd, acceptance.th))
  }
