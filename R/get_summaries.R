# Nothing here is exported


# summaries_omega ---------------------------------------------------------
summaries_omega <- function(bsbhaz_omega) {
  individuals <- nrow(bsbhaz_omega)
  iter <- ncol(bsbhaz_omega)
  means <- vector(mode = "double", length = individuals)
  prob_low <- vector(mode = "double", length = individuals)
  prob_high <- vector(mode = "double", length = individuals)
  acc_rate <- vector(mode = "double", length = individuals)
  for (i in 1:individuals) {
    means[[i]] <- mean(bsbhaz_omega[i, ])
    probs <- stats::quantile(bsbhaz_omega[i, ], probs = c(.025, .975))
    prob_low[[i]] <- probs[[1]]
    prob_high[[i]] <- probs[[2]]
    acc_rate[[i]] <- acceptance_rate(bsbhaz_omega[i, ])
  }
  data.frame("Individual" = 1:individuals,
             "Mean" = means,
             "Prob. Low 95%" = prob_low,
             "Prob. High 95%" = prob_high,
             "Acceptance Rate" = acc_rate,
             check.names = FALSE)
}


# summaries_lambda --------------------------------------------------------
summaries_lambda <- function(bsbhaz_lambda) {
  intervals <- nrow(bsbhaz_lambda)
  iter <- ncol(bsbhaz_lambda)
  means <- vector(mode = "double", length = intervals)
  prob_low <- vector(mode = "double", length = intervals)
  prob_high <- vector(mode = "double", length = intervals)
  for (i in 1:intervals) {
    means[[i]] <- mean(bsbhaz_lambda[i, ])
    probs <- stats::quantile(bsbhaz_lambda[i, ], probs = c(.025, .975))
    prob_low[[i]] <- probs[[1]]
    prob_high[[i]] <- probs[[2]]
  }
  data.frame("Interval" = 1:intervals,
             "Mean" = means,
             "Prob. Low 95%" = prob_low,
             "Prob. High 95%" = prob_high,
             check.names = FALSE)
}


# summaries_gamma ---------------------------------------------------------
summaries_gamma <- function(bsbhaz_gamma) {
  iter <- ncol(bsbhaz_gamma)
  means <- mean(bsbhaz_gamma[1, ])
  probs <- stats::quantile(bsbhaz_gamma[1, ], probs = c(.025, .975))
  prob_low <- probs[[1]]
  prob_high <- probs[[2]]
  acc_rate <- acceptance_rate(bsbhaz_gamma[1, ])
  data.frame("Gamma" = "Gamma",
             "Mean" = means,
             "Prob. Low 95%" = prob_low,
             "Prob. High 95%" = prob_high,
             "Acceptance Rate" = acc_rate,
             check.names = FALSE)
}


# summaries_theta ---------------------------------------------------------
summaries_theta <- function(bsbhaz_theta) {
  predictors <- nrow(bsbhaz_theta)
  iter <- ncol(bsbhaz_theta)
  means <- vector(mode = "double", length = predictors)
  prob_low <- vector(mode = "double", length = predictors)
  prob_high <- vector(mode = "double", length = predictors)
  acc_rate <- vector(mode = "double", length = predictors)
  for (i in 1:predictors) {
    means[[i]] <- mean(bsbhaz_theta[i, ])
    probs <- stats::quantile(bsbhaz_theta[i, ], probs = c(.025, .975))
    prob_low[[i]] <- probs[[1]]
    prob_high[[i]] <- probs[[2]]
    acc_rate[[i]] <- acceptance_rate(bsbhaz_theta[i, ])
  }
  data.frame("Predictor" = rownames(bsbhaz_theta),
             "Coefficient Mean" = means,
             "Prob. Low 95%" = prob_low,
             "Prob. High 95%" = prob_high,
             "Acceptance Rate" = acc_rate,
             check.names = FALSE)
}


# summaries_surv ----------------------------------------------------------

summaries_surv <- function(bsbhaz_surv) {
  times <- nrow(bsbhaz_surv)
  iter <- ncol(bsbhaz_surv)
  means <- vector(mode = "double", length = times)
  prob_low <- vector(mode = "double", length = times)
  prob_high <- vector(mode = "double", length = times)
  for (i in 1:times) {
    means[[i]] <- mean(bsbhaz_surv[i, ])
    probs <- stats::quantile(bsbhaz_surv[i, ], probs = c(.025, .975))
    prob_low[[i]] <- probs[[1]]
    prob_high[[i]] <- probs[[2]]
  }
  data.frame("t" = rownames(bsbhaz_surv),
             "S(t)" = means,
             "Prob. Low 95%" = prob_low,
             "Prob. High 95%" = prob_high,
             check.names = FALSE)
}
