# Nothing here is exported

# log_density_gamma -------------------------------------------------------

# Computes the proportional log-density of gamma
log_density_gamma <- function(gamma, proposal, omega1, omega2, y, gamma_d) {

  n <- length(omega1)
  sum_y <- sum(y)
  sum_omega <- sum(omega1) + sum(omega2)

  pi <- (2 * n - 2) * (log(1 + proposal) - log(1 + gamma)) +
    sum_y * (log(proposal) - log(gamma) + log(1 + proposal) - log(1 + gamma)) -
    sum_omega * (proposal - gamma)
  q <- log(proposal) - log(gamma)

  pi + q

}


# sample_gamma ------------------------------------------------------------

# Samples from gamma posterior
sample_gamma <- function(gamma, omega1, omega2, y, gamma_d) {

  if (is.null(gamma_d)) gamma_d <- 2
  proposal <- stats::rgamma(n = 1, shape = gamma_d, rate = gamma_d / gamma)
  l_density <- log_density_gamma(gamma, proposal, omega1, omega2, y, gamma_d)
  alpha <- min(exp(l_density), 1)
  u <- stats::runif(n = 1)

  gamma + (proposal - gamma) * (u <= alpha)

}
