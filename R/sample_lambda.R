# Nothing here is exported

# lambda_restriction ------------------------------------------------------

#Computes the restriction imposed by a single individual. For use in purrr::map.
lambda_restriction <- function(t,
                               omega,
                               x,
                               part_loc,
                               t_partition_low,
                               lambda_index,
                               theta,
                               lambda,
                               part_len) {

  omega_h <- omega * exp(-(x %*% theta))

  if (lambda_index < part_loc) {
    current_bound <- 0

    if (part_loc > 1) {
      for (k in 1:(part_loc - 1)) {
        current_bound <-
          current_bound + (lambda[k] * part_len) * (k != lambda_index)
      }
    }

    current_bound <-
      current_bound + lambda[part_loc] * (t - t_partition_low)
    current_bound <- (omega_h - current_bound) / part_len

    if (is.na(current_bound)) {
      warning("Lambda bound failed")
      current_bound <- 2
    }

  } else if (lambda_index == part_loc) {
    current_bound <- 0

    if (part_loc > 1) {
      for (k in 1:(part_loc - 1)) {
        current_bound <- current_bound + (lambda[k] * part_len)
      }
    }

    current_bound <- (omega_h - current_bound) / (t - t_partition_low)

    if (is.na(current_bound)) {
      warning("Lambda bound failed")
      current_bound <- 2
    }

  } else {

    current_bound <- 2

  }

  return(current_bound)

}


# get_min_bound -----------------------------------------------------------

# Computes the most restrictive bound for a given lambda
get_min_bound <- function(t,
                          omega,
                          x,
                          part_loc,
                          t_partition_low,
                          l_index,
                          theta,
                          lambda,
                          part_len) {

  bounds <-purrr::pmap_dbl(
    list(t, omega, x, part_loc, t_partition_low),
    function(x1, x2, x3, x4, x5) {
      lambda_restriction(x1, x2, x3, x4, x5, l_index, theta, lambda, part_len)
    }
  )

  return(min(bounds))

}


# sample_lambda -----------------------------------------------------------

# Samples a single scalar observation. For use in purrr::pmap.
sample_lambda <- function(u1, u2, alpha, beta, c1, c2, min_bound, part_count) {

  # u2 = 0 if lambda is the first interval
  # c2 = c1 if lambda is not in the first interval
  alpha_l <- alpha + u1 + u2 + part_count
  beta_l <- beta + c1 + c2
  
  if (min_bound < 1e-5) min_bound <- 1e-5
  denominator <- stats::pgamma(min_bound * beta_l, shape = alpha_l, rate = 1)
  if (denominator < 1e-5) denominator <- 1e-5
  unif <- stats::runif(n = 1)

  f <- function(x) {
    stats::pgamma(x, shape = alpha_l, rate = 1) / denominator - unif
  }
  
  if (f(min_bound * beta_l) <= 0) {
    solution <- min_bound * beta_l
  } else {
    solution <- stats::uniroot(f, lower = 0, upper = min_bound * beta_l)$root
  }

  
  lambda_prueba <- solution / beta_l
  if (lambda_prueba < 1e-5) lambda_prueba <- 1e-5
  
  return(lambda_prueba)

}
