CUpdTheta <-
  function(theta, m, lambda.r, times, delta, K, covar, tao, var.theta.str, var.theta.ini, acceptance) {
    p <- length(theta)
    theta.upd <- theta
    m.upd <- m
    for(s in 1:p){
      theta.str <- theta.upd
      theta.str[s] <- rnorm(1, mean = theta.upd[s], sd = sqrt(var.theta.str))
      m.str <- CGaM(times, tao, K, covar, theta.str)
      pr <- dnorm(theta.str[s], mean = 0, sd = sqrt(var.theta.ini), log = T) -
        dnorm(theta.upd[s], mean = 0, sd = sqrt(var.theta.ini), log = T) +
        (theta.str[s] - theta.upd[s]) * sum(covar[delta==1,s]) +
        sum(lambda.r * (m.upd - m.str))
      if(log(runif(1)) <= pr) {
        theta.upd[s] <- theta.str[s]
        m.upd <- m.str
        acceptance[s] <- acceptance[s] + 1
      }
    }
    return(list(theta.upd, acceptance))
  }
