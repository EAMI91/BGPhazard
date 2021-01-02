# Nothing here is exported

# log_density_u -----------------------------------------------------------

# Computes the proportional log-density of u
log_density_u <- function(u, l_k, l_k1, alpha, beta, c, index) {
  # index = 0 if u is the last, 1 in any other case

  u * (log(c) + index * log(c + beta) + log(l_k) + index * log(l_k1)) -
    lgamma(u + 1) - lgamma(alpha + u) * index

}


# prob_u ------------------------------------------------------------------

# Creates the distribution function for u
prob_u <- function(l, l1, alpha, beta, c, index) {

  u <- vector(mode = "numeric", length = 5e2L)
  acum <- 0
  j <- 0
  prob <- 1

  while (j <= 5e2 & prob > 1e-6) {
    pi_j <- exp(log_density_u(j, l, l1, alpha, beta, c, index))
    if (is.nan(pi_j)) {
      break
      }
    prueba_acum <- acum + pi_j
    if (is.infinite(prueba_acum)) {
      break
      }
    if (prueba_acum == 0) {
      break
      }
    acum <- prueba_acum
    prob <- pi_j / acum
    u[j + 1] <- pi_j
    j <- j + 1
  }

  u <- u[1:j] / acum
  prob_fun <- cumsum(u)
  
  return(prob_fun)

}


# sample_u ----------------------------------------------------------------

# Samples a single scalar observation. For use in purrr::pmap.
sample_u <- function(l, l1, alpha, beta, c, index_indicator) {

  distribution <- prob_u(l, l1, alpha, beta, c, index_indicator)

  if(is.na(distribution[1])) {
    distribution <- c(1)
    warning("No distribution for U")
  }

  u <- stats::runif(n = 1)
  index <- 1

  while (u > distribution[index]) {
    index <- index + 1
  }

  y <- index - 1
  
  return(y)
  
}
