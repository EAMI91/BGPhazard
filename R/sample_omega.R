# Not exported

# Samples a single scalar observation. For use in purrr::pmap.
sample_omega <- function(omega, y, cum_h, x, theta, gamma, omega_d = NULL) {

  if (is.null(omega_d)) omega_d <- y + 1

  bound <- cum_h * exp(x %*% theta)
  l1 <- max(bound, omega - omega_d)
  l2 <- omega + omega_d
  proposal <- stats::runif(n = 1, min = min(l1, l2), max = max(l1, l2))

  if (omega <= bound) {
    omega_out <- proposal
    return(omega_out)
  }

  log_alpha <-
    y * (log(proposal) - log(omega)) - (1 + gamma) * (proposal - omega)
  prob <- min(exp(log_alpha), 1)
  u <- stats::runif(n = 1)
  omega_out <- omega + (proposal - omega) * (u <= prob)
  
  return(omega_out)

}
