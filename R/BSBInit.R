#' Initial setup for BSBHaz model
#'
#' \code{BSBInit} creates the necessary data structure for use in
#' \code{\link{BSBHaz}}.
#'
#' This function reads and formats censored bivariate survival data in the
#' following way. If \code{df} is provided, failure times and censoring
#' indicadors are assumed to be columns named 't1', 't2', 'delta1', and
#' 'delta2'. Other columns not named 'id' (ignoring case) are taken to be
#' predictors. If \code{df} has no columns 'delta1' or 'delta2', observations
#' are taken as exact.
#'
#' If \code{df} is not provided, then \code{t1} and \code{t2} are expected to be
#' objects of class 'Surv' created by \code{\link[survival]{Surv}} and the model
#' does not use predictors. Only right-censored observations are supported. Only
#' \code{df} or \code{t1} and \code{t2} must be supplied. \code{df} argument
#' comes first for use in pipes.
#'
#' @param df A data frame with columns 't1', 't2', 'delta1', 'delta2'. Any other
#'   columns not named 'id' are taken to be predictors. These predictors must be
#'   numeric, i.e., \strong{categorical predictors must be one-hot encoded}.
#' @param t1,t2 Objects of class 'Surv' as created by
#'   \code{\link[survival]{Surv}}.
#' @param alpha,beta,c Doubles. Parameters for Markov gamma hazard priors.
#' @param part_len A double that gives the length of time partition intervals.
#' @param seed Random seed for variable initialization.
#'
#' @return An object of class '\code{BSBinit}'
#'
#' @examples
#' t1 <- survival::Surv(c(1, 2, 3))
#' t2 <- survival::Surv(c(1, 2, 3))
#'
#' init <- BSBInit(t1 = t1, t2 = t2, seed = 0)
#'
#' @export
BSBInit <- function(df = NULL,
                    t1 = NULL,
                    t2 = NULL,
                    alpha = 0.001,
                    beta = 0.001,
                    c = 1000,
                    part_len = 1,
                    seed = 42) {
  
  if (!is.null(df)) {
    if (!inherits(df, "data.frame")) {
      stop("df must be of class 'data.frame'")
    }
    if (!(is.null(t1) & is.null(t2))) {
      stop("Only one of df or t1 and t2 must be supplied")
    }
  }
  
  if (!is.null(t1)) {
    if (!inherits(t1, "Surv")) {
      stop("t1 is not of class 'Surv'")
    }
    if (attributes(t1)$type != "right") {
      stop(paste("t1 has censoring of type:", attributes(t1)$type))
    }
    stopifnot(attributes(t1)$type == "right")
    if (!is.null(df)) {
      stop("Only one of df or t1 and t2 must be supplied")
    }
    if (is.null(t2)) {
      stop("Missing t2")
    }
    if (length(t1) != length(t2)) {
      stop("t1 and t2 must have the same length")
    }
  }
  
  if (!is.null(t2)) {
    if (!inherits(t2, "Surv")) {
      stop("t2 is not of class 'Surv'")
    }
    if (attributes(t2)$type != "right") {
      stop(paste("t2 has censoring of type:", attributes(t2)$type))
    }
    if (!is.null(df)) {
      stop("Only one of df or t1 and t2 must be supplied")
    }
    if (is.null(t1)) {
      stop("Missing t1")
    }
    if (length(t1) != length(t2)) {
      stop("t1 and t2 must have the same length")
    }
  }
  
  # Handling the 'Surv' objects
  if (!is.null(t1)) {
    n_obs <- attributes(t1)$dim[[1]]
    delta1 <- as.double(t1)[(n_obs + 1):(n_obs * 2)]
    t1 <- as.double(t1)[1:n_obs]
    delta2 <- as.double(t2)[(n_obs + 1):(n_obs * 2)]
    t2 <- as.double(t2)[1:n_obs]
    has_predictors <- FALSE
    pred_matrix <- matrix(rep(0, times = n_obs), nrow = n_obs)
  }
  
  if (!is.null(df)) {
    stopifnot(!is.null(df[["t1"]]))
    stopifnot(!is.null(df[["t2"]]))
    df[["id"]] <- NULL
    df[["Id"]] <- NULL
    df[["iD"]] <- NULL
    df[["ID"]] <- NULL
    n_obs <- length(df[["t1"]])
    if (is.null(df[["delta1"]])) {
      delta1 <- rep(1, times = n_obs)
    } else {
      stopifnot(is.double(df[["delta1"]]) | is.integer(df[["delta1"]]))
      delta1 <- as.double(df[["delta1"]])
    }
    if (is.null(df[["delta2"]])) {
      delta2 <- rep(1, times = n_obs)
    } else {
      stopifnot(is.double(df[["delta2"]]) | is.integer(df[["delta2"]]))
      delta2 <- as.double(df[["delta2"]])
    }
    stopifnot(is.double(df[["t1"]]) | is.integer(df[["t1"]]))
    t1 <- as.double(df[["t1"]])
    stopifnot(is.double(df[["t2"]]) | is.integer(df[["t2"]]))
    t2 <- as.double(df[["t2"]])
    df[["t1"]] <- NULL
    df[["t2"]] <- NULL
    df[["delta1"]] <- NULL
    df[["delta2"]] <- NULL
    
    # Predictors
    if (length(names(df)) != 0) {
      has_predictors <- TRUE
      pred_matrix <- as.matrix(df)
      colnames(pred_matrix) <- colnames(df)
    } else {
      has_predictors <- FALSE
      pred_matrix <- matrix(rep(0, times = n_obs), nrow = n_obs)
    }
  }
  
  # Time partition
  max_t <- max(max(t1), max(t2))
  t_part <- partition(t = max_t, int_len = part_len)
  
  # Variable initialization
  set.seed(seed)
  rho <- stats::runif(n = 1)
  gamma <- rho / (1 - rho)
  omega1 <- stats::rgamma(n = n_obs, shape = 2, rate = 1)
  y <- stats::rpois(n = n_obs, lambda = omega1)
  omega2 <- stats::rgamma(n = n_obs, shape = 2, rate = 1)
  theta <- stats::rnorm(n = ncol(pred_matrix))
  n_intervals <- length(t_part) - 1
  lambda1 <- stats::rgamma(n = n_intervals, shape = alpha, rate = beta)
  lambda2 <- stats::rgamma(n = n_intervals, shape = alpha, rate = beta)
  lambda1 <- pmax(lambda1, rep(1e-5, times = n_intervals))
  lambda2 <- pmax(lambda2, rep(1e-5, times = n_intervals))
  u1 <- stats::rpois(n = n_intervals, lambda = lambda1)
  u2 <- stats::rpois(n = n_intervals, lambda = lambda2)
  
  list_out <- list(
    "t1" = t1, "t2" = t2, "delta1" = delta1, "delta2" = delta2,
    "alpha" = alpha, "beta" = beta, "c" = c,
    "gamma" = gamma, "omega1" = omega1, "omega2" = omega2, "y" = y,
    "theta" = theta, "pred_matrix" = pred_matrix, "u1" = u1, "u2" = u2,
    "lambda1" = lambda1, "lambda2" = lambda2, "t_part" = t_part
  )
  new_BSBinit(
    l = list_out,
    individuals = as.integer(n_obs),
    intervals = as.integer(n_intervals),
    has_predictors = has_predictors
  )
  
}
