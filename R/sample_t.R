# Nothing here is exported

# sample_t ----------------------------------------------------------------

# Samples censored observations from distribution function. For use in pmap.
sample_t <- function(t_orig,
                     t_prev,
                     omega,
                     delta,
                     max_part,
                     x,
                     theta,
                     partition,
                     lambda) {

  bound <- t_orig

  if (delta == 1) {
    return(t_orig)
  }

  u <- stats::runif(n = 1)

  f <- function(var) {
    cum_h(var, partition, lambda) * exp(x %*% theta) - u * omega
  }

  up <- max_part

  while (f(up) < 0) {
    up <- up + 10
  }

  proposal <- stats::uniroot(f, lower = 0 , upper = up)$root
  out <-
    t_prev + (proposal - t_prev) * (proposal > bound & proposal <= max_part)
  
  return(out)

}
