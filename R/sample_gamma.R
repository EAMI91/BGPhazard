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
  # q <- log(proposal) - log(gamma)

  return(pi)

}


# sample_gamma ------------------------------------------------------------

# Samples from gamma posterior
sample_gamma <- function(gamma, omega1, omega2, y, gamma_d) {

  # if (is.null(gamma_d)) gamma_d <- 2
  if (is.null(gamma_d)) gamma_d <- gamma * 0.5
  l1 <- max(0, gamma - gamma_d)
  l2 <- gamma + gamma_d
  # proposal <- stats::rgamma(n = 1, shape = gamma_d, rate = gamma_d / gamma)
  proposal <- stats::runif(n = 1, min = min(l1, l2), max = max(l1, l2))
  l_density <- log_density_gamma(gamma, proposal, omega1, omega2, y, gamma_d)
  alpha <- min(exp(l_density), 1)
  u <- stats::runif(n = 1)

  out <- gamma + (proposal - gamma) * (u <= alpha)
  
  return(out)

}
