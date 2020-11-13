# Nothing here is exported

# theta_restriction -------------------------------------------------------

# Computes the restriction imposed on a single theta by a single individual
theta_restriction <- function(t, omega, cum_h, x, theta_index, theta) {

  current_th <- theta[theta_index]
  current_x <- x[theta_index]

  if (current_x == 0) {
    bound <- 1e4
    return(bound)
  }

  out <-
    (log(omega) - log(cum_h) - x %*% theta + current_th * current_x) / current_x
  
  return(out)

  }


# get_min_bound_theta -----------------------------------------------------

# Computes the most restrictive bound for a single theta
get_min_bound_theta <- function(theta_index, t, omega, cum_h, x, theta) {

  bounds <- purrr::pmap_dbl(
    list(t, omega, cum_h, x),
    function(x1, x2, x3, x4) {
      theta_restriction(
        x1,
        x2,
        x3,
        x4,
        theta_index,
        theta
      )
    }
  )

  return(min(bounds))

}


# sample_theta ------------------------------------------------------------

# Samples a single posterior observation. For use in purrr::pmap.
sample_theta <- function(bound, sum_x, theta, theta_d = NULL) {

  if (is.null(theta_d)) theta_d <- 0.5 * theta
  l1 <- theta - theta_d
  l2 <- min(bound, theta + theta_d)
  proposal <- stats::runif(n = 1, min = min(l1, l2), max = max(l1, l2))

  if (theta > bound) {
    out <- proposal
    return(out)
  }

    l_rho <-
      proposal * (2 * sum_x - 0.5 * proposal) - theta * (2 * sum_x - 0.5 * theta)
    alpha <- min(exp(l_rho), 1)
    u <- stats::runif(n = 1)
    out <- theta + (proposal - theta) * (u <= alpha)
    
    return(out)

}
