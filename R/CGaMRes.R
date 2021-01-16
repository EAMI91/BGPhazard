#' Markov Gamma Model with Covariates
#' 
#' Posterior inference for the Bayesian non-parametric Markov gamma model with
#' covariates in survival analysis.
#' 
#' Computes the Gibbs sampler with the full conditional distributions of
#' Lambda and Theta (Nieto-Barajas, 2003) and arranges the resulting Markov
#' chain into a matrix which can be used to obtain posterior summaries. Prior
#' distributions for the re gression coefficients (Theta) are assumed independent normals
#' with zero mean and variance \code{var.theta.ini}.
#' 
#' @param data Double tibble. Contains failure times in the first column,
#' status indicator in the second, and, from the third to the last column, the
#' covariate(s).
#' @param type.t Integer. 1=computes uniformly-dense intervals; 2=length
#' intervals defined by user and 3=same length intervals.
#' @param length Integer. Interval length of the partition.
#' @param K Integer. Partition length for the hazard function.
#' @param alpha Nonnegative entry vector. Small entries are recommended in
#' order to specify a non-informative prior distribution.
#' @param beta Nonnegative entry vector. Small entries are recommended in order
#' to specify a non-informative prior distribution.
#' @param c.r Nonnegative vector. The higher the entries, the higher the correlation of 
#' two consecutive intervals.
#' @param c.nu Tuning parameter for the proposal distribution for c.
#' @param var.theta.str Double. Variance of the proposal normal distribution
#' for theta in the Metropolis-Hastings step.
#' @param var.theta.ini Double. Variance of the prior normal distribution for theta.
#' @param a.eps Double. Shape parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param b.eps Double. Scale parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param type.c 1=defines \code{c.r} as a zero-entry vector; 2=lets the user
#' define \code{c.r} freely; 3=assigns \code{c.r} by computing an exponential
#' prior distribution with mean 1; 4=assigns \code{c.r} an exponential hierarchical
#' distribution with mean \code{epsilon} which in turn has a Ga(a.eps, b.eps)
#' distribution.
#' @param epsilon Double. Mean of the exponential distribution assigned to
#' \code{c.r} when \code{type.c = 3}.
#' @param iterations Integer. Number of iterations including the \code{burn.in}
#' to be computed for the Markov chain.
#' @param burn.in Integer. Length of the burn-in period for the Markov chain.
#' @param thinning Integer. Factor by which the chain will be thinned. Thinning
#' the Markov chain reduces autocorrelation.
#' @param printtime Logical. If \code{TRUE}, prints out the execution time.
#' @note It is recommended to verify chain's stationarity. This can be done by
#' checking each element individually. See \link{CGaPlotDiag}
#' To obtain posterior summaries of the coefficients use function
#' \link{CGaPloth}.
#' @seealso \link{CGaPlotDiag}, \link{CGaPloth}
#' @references - Nieto-Barajas, L. E. (2003). Discrete time Markov gamma
#' processes and time dependent covariates in survival analysis. \emph{Bulletin
#' of the International Statistical Institute 54th Session}. Berlin. (CD-ROM).
#' 
#' - Nieto-Barajas, L. E. & Walker, S. G. (2002). Markov beta and gamma
#' processes for modelling hazard rates. \emph{Scandinavian Journal of
#' Statistics} \strong{29}: 413-424.
#' @examples
#' 
#' 
#' 
#' ## Simulations may be time intensive. Be patient.
#' 
#' ## Example 1
#' #  data(leukemiaFZ)
#' #  leukemia1 <- leukemiaFZ
#' #  leukemia1$wbc <- log(leukemiaFZ$wbc)
#' #  CGEX1 <- CGaMRes(data = leukemia1, K = 10, iterations = 100, thinning = 1)
#' 
#' ## Example 2. Refer to "Cox-gamma model example" section in package vignette for details.
#' # SampWeibull <- function(n, a = 10, b = 1, beta = c(1, 1)) {
#' #   M <- tibble(i = seq(n), x_i1 = runif(n), x_i2 = runif(n), 
#' #               t_i = rweibull(n, shape = b, 
#' #                                 scale = 1 / (a * exp(x_i1*beta[1] + x_i2*beta[2]))),
#' #               c_i = rexp(n), delta = t_i > c_i,
#' #               `min{c_i, d_i}` = min(t_i, c_i))
#' #   return(M)
#' # }
#' #  dat <- SampWeibull(100, 0.1, 1, c(1, 1))
#' #  dat <- dat %>% select(4,6,2,3) 
#' #  CG <- CGaMRes(data = leukemia1, K = 10, iterations = 100, thinning = 1)
#' #  CGaPloth(CG)
#' 
#' 
#' 
#' @export CGaMRes
CGaMRes <-
  function(data, type.t = 2, length = 1, K = 5, alpha = rep(0.01, K),
           beta = rep(0.01, K), c.r = rep(1, K - 1), c.nu = 1, 
           var.theta.str = 25, var.theta.ini = 100,
           a.eps = 0.1, b.eps = 0.1, 
           type.c = 4, epsilon = 1, iterations = 1000, 
           burn.in = floor(iterations * 0.2), thinning = 3, printtime = TRUE) {
    tInit <- proc.time()
    data <- tibble::as_tibble(data)
    print(sprintf("Using %s as times and %s as delta, status indicator. The other variables are used as covariables",names(data)[1], names(data)[2]))
    times <- as.numeric(dplyr::pull(data, 1))
    delta <- as.numeric(dplyr::pull(data, 2))
    covar <- dplyr::select(data, -c(1, 2))
    
    covar2 <- covar
    
    median.obs <- purrr::map_df(covar, ~quantile(x = .x, probs = .5))
    
    k.const <- as.numeric(covar %>% dplyr::summarise_all(.funs = ~max(abs(.x))))
    covar %<>% purrr::modify2(.y = k.const, .f = ~.x/.y)
    
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
      if(is.null(utao)) stop("If type.t = 2 you need to specify utao.")
      utao <- sort(utao)
      if(utao[1]!=0){
        warning("The first value of the partition needs to be 0, utao fixed and now starting with 0.")
        utao <- c(0, utao)
      } 
      if(max(times) > max(utao)){
        utao <- c(utao,max(times))
        warning("The last value of the partition needs to be", max(times),", utao fixed and set to ",max(times),".")
      }
      K <- length(utao) - 1
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
        stop (c("Invalid argument: 'c.r' must have length, "), K - 1)
      }
      if (sum(abs(c.r - round(c.r)) > tol) != 0 || min(c.r) < 0) {
        stop ("Invalid argument: 'c.r' entries must be nonnegative integers.")
      }
    }
    if (type.c == 1 && sum(abs(c.r)) != 0 ) {
      c.r <- rep(0, K - 1)
      warning (c("'c.r' redefined as rep(0,", K - 1, ") because type.c = 1."))
    }
    if ((type.c == 3 || type.c == 4) && epsilon < 0) {
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
    if (class(thinning) != "numeric") {
      stop ("Invalid argument: 'thinning' must be a numeric value.")
    }
    if (thinning <= 0 || abs(thinning - round(thinning)) > tol 
        || thinning > 0.1 * iterations) {
      stop ("Invalid argument: 'thinning' must be a postitive integer smaller than 
            iterations * 0.10 = ", iterations * 0.1, ".")
    }
    if (printtime != TRUE && printtime != FALSE) {
      stop ("Invalid argument: 'printtime' must be a logical value.")
    }
    tao <- Tao(times, delta, type.t, K, length)
    t.unc <- sort(times[delta==1])
    n <- readr::parse_integer(as.character(table(cut(t.unc,tao))))
    if (type.c == 3) {
      c.r <- rep(5, (K - 1))
    }
    if (type.c == 4) {
      Epsilon <- rep(NA, iterations)
    }  
    p <- ncol(covar)
    acceptance.th <- rep(0,p)
    acceptance.c <- 0
    Theta <- matrix(NA, nrow = iterations, ncol = p)
    Lambda <- matrix(NA, nrow = iterations, ncol = K)
    U <- matrix(NA, nrow = iterations, ncol = K - 1)
    C <- matrix(NA, nrow = iterations, ncol = K - 1)
    lambda.r <- rep(0.1, K)
    theta <- rep(1, p)
    cat(paste("Iterating...", "\n"), sep = "")
    pb <- dplyr::progress_estimated(iterations)
    covar <- as.matrix(covar)
    for(j in seq_len(iterations)) {
      pb$tick()$print()
      u.r <- UpdU(alpha, beta, c.r, lambda.r)
      m <- CGaM(times, tao, K, covar, theta)
      lambda.r <- UpdLambda(alpha, beta, c.r, u.r, n, m)
      aux.th <- CUpdTheta(theta, m, lambda.r, times, delta, K, covar, tao, var.theta.str, var.theta.ini, acceptance.th)
      theta <- aux.th[[1]]
      acceptance.th <- aux.th[[2]]
      if (type.c == 3 || type.c == 4) {
        if (type.c == 4) {
          epsilon <- rgamma(1, shape = a.eps + K, scale = 1 / (b.eps + sum(c.r)))
        }
        auxc.r <- GaUpdC(alpha, beta, c.r, lambda.r, u.r, epsilon, c.nu, acceptance.c)
        c.r <- auxc.r[[1]]
        acceptance.c <- auxc.r[[2]]
      }
      Lambda[j, ] <- lambda.r
      U[j, ] <- u.r
      C[j, ] <- c.r
      if (type.c == 4) Epsilon[j] <- epsilon
      Theta[j, ] <- theta
    }
    Lambda <- Lambda[seq(burn.in + 1, iterations, thinning), ]
    U <- U[seq(burn.in + 1, iterations, thinning), ]
    C <- C[seq(burn.in + 1, iterations, thinning), ]
    if (type.c == 4) Epsilon <- Epsilon[seq(burn.in + 1, iterations, thinning)]
    Theta <- Theta[seq(burn.in + 1, iterations, thinning), ]
    Theta <- sweep(Theta, MARGIN=2,k.const, `/`)
    
    
    
    rows <- nrow(Lambda)
    aux.median.obs <- as.numeric(median.obs)
    
    Lambda.median.obs <- tibble::as_tibble(Lambda)
    eff <- as.numeric(exp(Theta%*%aux.median.obs))
    Lambda.median.obs <- dplyr::mutate_all(Lambda.median.obs,.f = ~.x*eff)
    
    ss <- max(tao) * seq.int(0,100) / 100
    X <- as.matrix(unname(Lambda))
    
    writeLines(c("","Done.","Generating baseline survival function estimates."))
    pb <- dplyr::progress_estimated(length(ss))
    S <-  purrr::map(ss, function(s = .x){
      pb$tick()$print()
      do.call(base::c,purrr::map(seq_len(rows),.f= ~exp(-sum((s > tao[-1]) * tao[-1] * X[.x,] + 
                                               (s > tao[-length(tao)] & s <= tao[-1]) * s * X[.x,] -
                                               (s > tao[-length(tao)]) * tao[-(length(tao))] * X[.x,])
      )))
    })
    
    X.median.obs <- as.matrix(unname(Lambda.median.obs))
    writeLines(c("","Done.","Generating survival function estimates of the median observation."))
    pb <- dplyr::progress_estimated(length(ss))
    S.median.obs <-  purrr::map(ss, function(s = .x){
      pb$tick()$print()
      do.call(base::c,purrr::map(seq_len(rows),.f= ~exp(-sum((s > tao[-1]) * tao[-1] * X.median.obs[.x,] + 
                                                        (s > tao[-length(tao)] & s <= tao[-1]) * s * X.median.obs[.x,] -
                                                        (s > tao[-length(tao)]) * tao[-(length(tao))] * X.median.obs[.x,])
      )))
    })
    # H[2:rows, 2:101] <- -log(S[2:rows, 2:101])
    S <- purrr::map(.x = 1, ~S)
    S.median.obs <- purrr::map(.x = 1, ~S.median.obs)
    Lambda <- purrr::map(.x = 1, ~tibble::as_tibble(Lambda))
    Lambda.median.obs <- purrr::map(.x = 1, ~Lambda.median.obs) 
    cat(paste("Done.", "\n"), sep = "")
    if (printtime) {
      cat(">>> Total processing time (sec.):\n")
      print(procTime <- proc.time() - tInit)
    }
    if(type.c == 4) {
      X = tibble::enframe(list(Lambda = Lambda, Lambda.m = Lambda.median.obs, 
                       U = U, C = C, Epsilon = Epsilon, Theta = Theta)) } else { 
                         X = tibble::enframe(list(Lambda = Lambda, Lambda.m = Lambda.median.obs, U = U, C = C, Theta = Theta))
                       }
    out <- tibble::enframe(list(times = times, delta = delta, data = covar2, type.t = type.t, 
                        tao = tao, K = K, t.unc = t.unc, iterations = rows, burn.in = burn.in, thinning = thinning, 
                        acceptance = tibble::enframe(list(a.th = acceptance.th/iterations, a.c = acceptance.c/((K-1)*iterations))),
                        simulations = X, p = p, s = ss, S = S, S.m = S.median.obs
                        ))
    return(out)
  }
