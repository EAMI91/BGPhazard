#' Bayesian Semiparametric Cure Rate Model with an Unknown Threshold and
#' Covariate Information
#' 
#' Posterior inference for the bayesian semiparmetric cure rate model with
#' covariates in survival analysis.
#' 
#' Computes the Gibbs sampler with the full conditional distributions of
#' all model parameters (Nieto-Barajas & Yin, 2008) and arranges the resulting Markov
#' chain into a tibble which can be used to obtain posterior summaries. Prior
#' distributions for the regression coefficients Theta and Delta are assumend
#' independent normals with zero mean and variance \code{var.theta.ini},
#' \code{var.delta.ini}, respectively.
#' 
#' @param data Double tibble. Contains failure times in the first column,
#' status indicator in the second, and, from the third to the last column, the
#' covariate(s).
#' @param covs.x Character. Names of covariables to be part of the
#' multiplicative part of the hazard
#' @param covs.y Character. Names of covariables to determine the cure
#' threshold por each patient.
#' @param type.t Integer. 1=computes uniformly-dense intervals; 2=length
#' intervals defined by the user and 3=same length intervals.
#' @param length Integer. Interval length of the partition.
#' @param K Integer. Partition length for the hazard function.
#' @param alpha Nonnegative entry vector. Small entries are recommended in
#' order to specify a non-informative prior distribution.
#' @param beta Nonnegative entry vector. Small entries are recommended in order
#' to specify a non-informative prior distribution.
#' @param c.r Nonnegative vector. The higher the entries, the higher the correlation of two consective intervals.
#' @param c.nu Tuning parameter for the proposal distribution for c. 
#' Only when \code{type.c} is 3 or 4.
#' @param var.theta.str Double. Variance of the proposal normal distribution
#' for theta in the Metropolis-Hastings step.
#' @param var.delta.str Double. Variance of the proposal normal distribution
#' for delta in the Metropolis-Hastings step.
#' @param var.theta.ini Double. Variance of the prior normal distribution for theta.
#' @param var.delta.ini Double. Variance of the prior normal distribution for delta.
#' from the acceptance ratio in the Metropolis-Hastings algorithm for delta*.
#' @param type.c 1=defines \code{c.r} as a zero-entry vector; 2=lets the user
#' define \code{c.r} freely; 3=assigns \code{c.r} an exponential prior
#' distribution with mean 1; 4=assigns \code{c.r} an exponential hierarchical 
#' distribution with mean \code{epsilon} which in turn has a Ga(a.eps, b.eps)
#' distribution.
#' @param a.eps Double. Shape parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param b.eps Double. Scale parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param epsilon Double. Mean of the exponencial distribution assigned to
#' \code{c.r} when \code{type.c = 3}.
#' @param iterations Integer. Number of iterations including the \code{burn.in}
#' to be computed for the Markov chain.
#' @param burn.in Integer. Length of the burn-in period for the Markov chain.
#' @param thinning Integer. Factor by which the chain will be thinned. Thinning
#' the Markov chain reduces autocorrelation.
#' @param printtime Logical. If \code{TRUE}, prints out the execution time.
#' @note It is recommended to verify chain's stationarity. This can be done by
#' checking each element individually. See \code{\link{CCuPlotDiag}}.
#' @seealso \link{CCuPlotDiag}, \link{CCuPloth}
#' @references - Nieto-Barajas, L. E., & Yin, G. (2008). Bayesian
#' semiparametric cure rate model with an unknown threshold. Scandinavian
#' Journal of Statistics, 35(3), 540-556.
#' https://doi.org/10.1111/j.1467-9469.2007.00589.x
#' 
#' - Nieto-barajas, L. E. (2002). Discrete time Markov gamma processes and time
#' dependent covariates in survival analysis. Statistics, 2-5.
#' @examples
#' 
#' 
#' 
#' # data(BMTKleinbook)
#'     # res <- CCuMRes(BMTKleinbook, covs.x = c("tTransplant","hodgkin","karnofsky","waiting"),
#'     #                covs.y = c("tTransplant","hodgkin","karnofsky","waiting"),
#'     #                        type.t = 2, K = 72, length = 30,
#'     #                        alpha = rep(2,72), beta = rep(2,72), c.r = rep(50, 71), type.c = 2,
#'     #                        var.delta.str = .1, var.theta.str = 1,
#'     #                        var.delta.ini = 100, var.theta.ini = 100,
#'     #                        iterations = 100, burn.in = 10, thinning = 1)
#' 
#' 
#' 
#' @export CCuMRes
CCuMRes <-
  function(data, covs.x = names(data)[seq.int(3,ncol(data))], 
           covs.y = names(data)[seq.int(3,ncol(data))], 
           type.t = 3, length, K = 50, alpha = rep(0.01, K),
           beta = rep(0.01, K), c.r = rep(0, K - 1), c.nu = 1, 
           var.theta.str = 25, var.delta.str = 25, var.theta.ini = 100, var.delta.ini = 100,
           type.c = 4, a.eps = 0.1, b.eps = 0.1, epsilon = 1, iterations = 5000, 
           burn.in = floor(iterations * 0.2), thinning = 3, printtime = TRUE) {
    tInit <- proc.time()
    data <- tibble::as_tibble(data)
    writeLines(c(sprintf("Using %s as times and %s as delta, status indicator.",names(data)[1], names(data)[2]),
                 "The other variables are used as covariables"))
    times <- as.numeric(dplyr::pull(data, 1))
    delta <- as.numeric(dplyr::pull(data, 2))
    covar <- dplyr::select(data, -c(1, 2))
    
    covar2 <- covar
    
    median.obs <- purrr::map_df(covar, ~quantile(x = .x, probs = .5))
    
    k.const <- as.numeric(covar %>% dplyr::summarise_all(.funs = ~max(abs(.x))))
    covar %<>% purrr::modify2(.y = k.const, .f = ~.x/.y)
    
    
    median.obs.x <- dplyr::select(median.obs, !!covs.x)
    median.obs.y <- dplyr::select(median.obs, !!covs.y)
    covs.x <- as.matrix(dplyr::select(covar, !!covs.x))
    covs.y <- as.matrix(dplyr::select(covar, !!covs.y))
    covar <- as.matrix(covar)
    
    if (min(times) < 0) {
      stop ("Invalid argument: 'times' must be a nonnegative vector.")
    }
    if (min((delta  ==  0) + (delta  ==  1 )) == 0) {
      stop ("Invalid argument: 'delta' must have 0 - 1 entries.")
    }
    if (length(times) != length(delta)) {
      stop ("Invalid argument: 'times' and 'delta' must have same length.")
    }
    if (type.t == 2) {
      m <- ceiling(max(times))
      if (length > m) {
        stop (c("type.t = 2 requires length <=", m))
      }
      K <- ceiling(ceiling(max(times))/length)
    }
    if (type.t == 1 || type.t == 3) {
      if (class(try(K != 0, TRUE)) == "try-error") {
        K.aux <- 5
        warning ("'K' value not specified. 'K' fixed at ", K.aux, ".")
      } else {K.aux <- K}
      K <- K.aux
    }
    tol <- .Machine$double.eps ^ 0.5
    if (abs(type.t - round(type.t)) > tol || type.t < 1 || type.t > 3) {
      stop ("Invalid argument: 'type.t' must be an integer between 1 and 3.")
    }
    if (K <= 2 || abs(K - round(K)) > tol) {
      stop ("Invalid argument: 'K' must be an integer greater than 2.")
    }
    if (length(alpha) != K || length(beta) != K) {
      stop (c("Invalid argument: 'alpha', 'beta', must have length "), K)
    }
    if (min(c(alpha, beta)) < 0) {
      stop ("Invalid argument: 'alpha' and 'beta' must have nonnegative entries.")
    } 
    if (abs(type.c - round(type.c)) > tol || type.c < 1 || type.c > 4) {
      stop ("Invalid argument: 'type.c' must be an integer between 1 and 4.")
    }
    if (type.c == 1 || type.c == 2) {
      if (length(c.r) != (K - 1)) {
        stop (c("Invalid argument: 'c.r' must have length, ", K - 1))
      }
      if (sum(abs(c.r - round(c.r)) > tol) != 0 || min(c.r) < 0) {
        stop ("Invalid argument: 'c.r' entries must be nonnegative integers.")
      }
    }
    if (type.c == 1 && sum(abs(c.r)) != 0 ) {
      c.r <- rep(0, K - 1)
      warning (c("'c.r' redefined as rep0,", K - 1, ") because type.c = 1."))
    }
    if (type.c == 3 && epsilon < 0) {
      stop ("Invalid argument: 'epsilon' must be nonnegative.")
    }
    if (iterations <= 0 || abs(iterations - round(iterations)) > tol 
        || iterations < 50) {
      stop ("Invalid argument: 'iterations' must be an integer greater than 50.")
    }
    if (burn.in < 0 || abs(burn.in - round(burn.in)) > tol 
        || burn.in > iterations*0.9) {
      stop ("Invalid argument: 'burn.in' must be a postitive integer smaller than 
            iterations = ", iterations * 0.9, ".")
    }
    if (class(thinning)!= "numeric") {
      stop ("Invalid argument: 'thinning' must be a logical value.")
    }
    if (thinning <= 0 || abs(thinning - round(thinning)) > tol 
        || thinning > 0.1 * iterations) {
      stop ("Invalid argument: 'thpar' must be a postitive integer smaller than 
            iterations * 0.10 = ", iterations * 0.1, ".")
    }
    if (printtime != TRUE && printtime != FALSE) {
      stop ("Invalid argument: 'printtime' must be a logical value.")
    }
    
    tao <- Tao(times, delta, type.t, K, length)
    
    t.unc <- sort(times[delta == 1])
    n <- readr::parse_integer(as.character(table(cut(t.unc,tao))))
    acceptance.c <- 0
    if (type.c %in% c(3,4)) {
      c.r <- rep(5, (K - 1))
      Epsilon <- rep(NA, iterations)
    }
    p <- ncol(covs.x)
    p2 <- ncol(covs.y)
    acceptance.th <- rep(0,p)
    acceptance.d <- rep(0,p2)
    ind <- nrow(covs.x)
    Theta <- matrix(NA, nrow = iterations, ncol = p)
    Lambda <- matrix(NA, nrow = iterations, ncol = K)
    U <- matrix(NA, nrow = iterations, ncol = K - 1)
    C <- matrix(NA, nrow = iterations, ncol = K - 1)
    Z <- matrix(NA, nrow = iterations, ncol = ind)
    Delta <- matrix(NA, nrow = iterations, ncol = p)
    k_i <- rep(1, length(times))
    k_i[delta==1] <- as.numeric(cut(times[delta==1],tao,labels = seq_len(length(tao)-1)))
    z <- k_i
    lambda.r <- rep(0.1, K)
    theta <- rep(0, p)
    delta.r <- rep(0, p2)
    cat(paste("Iterating...", "\n"), sep = "")
    pb <- dplyr::progress_estimated(iterations)
    for(j in seq_len(iterations)) {
      pb$tick()$print()
      u.r <- UpdU(alpha, beta, c.r, lambda.r)
      W <- CCuW(theta, times, K, covs.x, tao, ind)
      IDW <- purrr::map2(W, .y = z, ~(seq_len(K) <= .y)*.x)
      m <- purrr::reduce(IDW, `+`)
      
      z <- CCuUpdZ(times,tao, lambda.r, W, delta.r, covs.y, k_i)
      lambda.r <- UpdLambda(alpha, beta, c.r, u.r, n, m)
      aux.th <- CCuUpdTheta(theta, lambda.r, times, delta, K, covs.x, tao, ind, z, var.theta.str, var.theta.ini, acceptance.th)
      theta <- aux.th[[1]]
      acceptance.th <- aux.th[[2]]
      aux.d <- CCuUpdDelta(delta.r, covs.y, z, var.delta.str, var.delta.ini, acceptance.d)
      delta.r <- aux.d[[1]]
      acceptance.d <- aux.d[[2]]
      if (type.c %in% c(3,4)) {
        if (type.c == 4) {
          epsilon <- rgamma(1, shape = a.eps + K, scale = 1 / (b.eps + sum(c.r)))
        }
        auxc.r <- GaUpdC(alpha, beta, c.r, lambda.r, u.r, epsilon, c.nu, acceptance.c)
        c.r <- auxc.r[[1]]
        acceptance.c <- auxc.r[[2]]
      }
      C[j, ] <- c.r
      Lambda[j, ] <- lambda.r
      U[j, ] <- u.r
      Z[j, ] <- z
      Theta[j, ] <- theta
      Delta[j, ] <- delta.r
      if (type.c == 4) Epsilon[j] <- epsilon
    }
    Lambda <-  Lambda[seq(burn.in + 1, iterations, thinning), ]
    U <- U[seq(burn.in + 1, iterations, thinning), ]
    C <- C[seq(burn.in + 1, iterations, thinning), ]
    Z <- Z[seq(burn.in + 1, iterations, thinning), ]
    Theta <- Theta[seq(burn.in + 1, iterations, thinning), ]
    Theta <- sweep(Theta, MARGIN=2,k.const, `/`)
    
    Delta <- Delta[seq(burn.in + 1, iterations, thinning), ]
    Delta <- sweep(Delta, MARGIN=2,k.const, `/`)
    if (type.c == 4) Epsilon <- Epsilon[seq(burn.in + 1, iterations, thinning)]
    
    rows <- nrow(Lambda)
    
    aux.median.obs.x <- as.numeric(median.obs.x)
    aux.median.obs.y <- as.numeric(median.obs.y)
    
    writeLines(c("","Done.","Generating predictive values for Z por the median observation."))
    
    
    z_median.obs <- purrr::map_int(purrr::map(seq_len(nrow(Delta)),
                                ~exp(purrr::reduce(purrr::map2(.x = aux.median.obs.y, .y = Delta[.x,], .f = ~.x*.y), `+`))),
                            ~{rpois(n = 1, lambda = .x)}) + 1
    z_median.obs[z_median.obs>K] <- K
    
    writeLines(c("", "Done.", "Generating predictive hazard rates for the median observation."))
    
    Lambda.median.obs <- purrr::map_dfc(seq_len(ncol(Lambda)), ~ ((.x <= z_median.obs)*Lambda[,.x]))
    
    X <- as.matrix(unname(Lambda.median.obs))
    
    
    
    writeLines(c("","Done.","Generating cure rate for the median observation."))
    pb <- dplyr::progress_estimated(rows)
    Pi.m <- do.call(base::c, purrr::map(seq_len(rows), 
                                 .f = ~ {
                                   pb$tick()$print()
                                   exp(-sum(exp(sum(Theta[.x,] * aux.median.obs.x)) * (tao[-1] - tao[-length(tao)]) * 
                                              Lambda.median.obs[.x,])
                                   )}))
    
    writeLines(c("","Done.","Generating survival function estimates of the median observation."))
    
    ss <- max(tao) * seq.int(0, 100) / 100
    
    pb <- dplyr::progress_estimated(length(ss))
    
    S <- purrr::map_dfc(ss, function(s = .x){
      pb$tick()$print()
      do.call(base::c, purrr::map(seq_len(rows), .f= ~exp(-sum((s > tao[-1]) * tao[-1] * X[.x,] * exp(sum(Theta[.x,] * aux.median.obs.x)) +
                                                          (s > tao[-length(tao)] & s <= tao[-1]) * s * X[.x,] * exp(sum(Theta[.x,] * aux.median.obs.x)) -
                                                          (s > tao[-length(tao)]) * tao[-(length(tao))] * X[.x,] * exp(sum(Theta[.x,] * aux.median.obs.x)))
      )))
    })
    
    S <- purrr::map(.x = 1, ~S)
    eff <- as.numeric(exp(Theta%*%aux.median.obs.x))
    Lambda.median.obs <- dplyr::mutate_all(Lambda.median.obs,.f = ~.x*eff)
    Lambda.median.obs <- purrr::map(.x = 1, ~Lambda.median.obs)
    cat(c("\n","Done.", "\n"), sep = "")
    if (printtime) {
      cat(">>> Total processing time (sec.):\n")
      print(procTime <- proc.time() - tInit)
    }
    if(type.c == 4) {
      X = list(Lambda = tibble::as_tibble(Lambda), Lambda.m = Lambda.median.obs,
               U = U, C = C, Theta = Theta, Delta = Delta, Z.m = z_median.obs, Pi.m = Pi.m, Epsilon = Epsilon)} else { 
                 X = list(Lambda = tibble::as_tibble(Lambda), Lambda.m = Lambda.median.obs, U = U, C = C, Theta = Theta, Delta = Delta, Z = Z, Z.m = z_median.obs, Pi.m = Pi.m)
               }
    
    X <- tibble::enframe(X)
    
    out <- tibble::enframe(list(times = times, delta = delta, data = covar2, covs.x = covs.x, covs.y = covs.y, type.t = type.t, 
                        tao = tao, K = K, t.unc = t.unc, iterations = rows, burn.in = burn.in, thinning = thinning, 
                        acceptance = tibble::enframe(list(a.d = acceptance.d/iterations, a.th = acceptance.th/iterations, a.c = acceptance.c/((K-1)*iterations))),
                        simulations = X, p = p, s = ss, S = S))
    
    return(out)
  }
