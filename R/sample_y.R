# Nothing here is exported

# log_density_y -----------------------------------------------------------

# Computes the proportional log-density of y
log_density_y <- function(y, omega1, omega2, gamma) {

  y * (log(gamma) + log(1 + gamma) + log(omega1) + log(omega2)) -
    lgamma(1 + y) - lgamma(2 + y)

}


# prob_y ------------------------------------------------------------------

# Creates the distribution function for y
prob_y <- function(omega1, omega2, gamma) {

  y <- vector(mode = "numeric", length = 5e2)
  acum <- 0
  j <- 0
  prob <- 1

  while (j <= 5e2 & prob > 1e-6) {
    gj <- exp(log_density_y(j, omega1, omega2, gamma))
    if (is.nan(gj)) {
      break
    }
    prueba_acum <- acum + gj
    if (is.infinite(prueba_acum)) {
      break
    }
    if (prueba_acum == 0) {
      break
    }
    acum <- prueba_acum
    prob <- gj / acum
    y[j + 1] <- gj
    j <- j + 1
  }

  y <- y[1:j] / acum
  prob_fun <- cumsum(y)
  
  return(prob_fun)

}


# sample_y ----------------------------------------------------------------

# Samples a single scalar observation. For use in purrr::pmap.
sample_y <- function(omega1, omega2, gamma) {

  distribution <- prob_y(omega1, omega2, gamma)

  if(is.na(distribution[1])) {
    distribution <- c(1)
    warning("No distribution for Y")
  }

  u <- stats::runif(n = 1)
  index <- 1

  while (u > distribution[index]) {
    index <- index + 1
  }

  y <- index - 1
  
  return(y)

}
