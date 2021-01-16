#' Bayesian Semiparametric Cure Rate Model with an Unknown Threshold
#' 
#' Posterior inference for the bayesian semiparametric cure rate model in
#' survival analysis.
#' 
#' 
#' Computes the Gibbs sampler with the full conditional distributions of
#' all model parameters (Nieto-Barajas & Yin 2008) and arranges the resulting Markov
#' chain into a tibble which can be used to obtain posterior summaries. 
#'
#' 
#' @param times Numeric positive vector. Failure times.
#' @param delta Logical vector. Status indicator. \code{TRUE} (1) indicates
#' exact lifetime is known, \code{FALSE} (0) indicates that the corresponding
#' failure time is right censored.
#' @param type.t Integer. 1=computes uniformly-dense intervals; 2=
#' partition arbitrarily defined by the user with parameter utao and 3=same length intervals.
#' @param K Integer. Partition length for the hazard function if
#' \code{type.t}=1 or \code{type.t}=3.
#' @param utao vector. Partition specified by the user when type.t = 2. The first value of 
#' the vector has to be 0 and the last one the maximum observed time, either censored or uncensored.
#' @param alpha Nonnegative entry vector. Small entries are recommended in
#' order to specify a non-informative prior distribution.
#' @param beta Nonnegative entry vector. Small entries are recommended in order
#' to specify a non-informative prior distribution.
#' @param c.r Nonnegative vector. The higher the entries, the higher the correlation of two consecutive intervals.
#' @param type.c 1=defines \code{c.r} as a zero-entry vector; 2=lets the user
#' define \code{c.r} freely; 3=assigns \code{c.r} by computing an exponential
#' prior distribution with mean epsilon; 4=assigns \code{c.r} by computing an exponential hierarchical
#' distribution with mean \code{epsilon} which in turn has a Ga(a.eps, b.eps)
#' distribution.
#' @param epsilon Double. Mean of the exponential distribution assigned to
#' \code{c.r} when \code{type.c = 3}. When \code{type.c = 4}, \code{epsilon} is
#' assigned a Ga(a.eps,b.eps) distribution.
#' @param c.nu Tuning parameter for the proposal distribution for c.
#' @param a.eps Numeric. Shape parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param b.eps Numeric. Scale parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param a.mu Numeric. Shape parameter for the prior gamma distribution of
#' mu
#' @param b.mu Numeric. Scale parameter for the prior gamma distribution of
#' mu
#' @param iterations Integer. Number of iterations including the \code{burn.in}
#' to be computed for the Markov Chain.
#' @param burn.in Integer. Length of the burn-in period for the Markov chain.
#' @param thinning Integer. Factor by which the chain will be thinned. Thinning
#' the Markov chain is to reduces autocorrelation.
#' @param printtime Logical. If \code{TRUE}, prints out the execution time.
#' @note It is recommended to verify chain's stationarity. This can be done by
#' checking each element individually. See \code{\link{CuPlotDiag}}.

#' @examples
#' 
#' 
#' ## Simulations may be time intensive. Be patient.
#' ## Example 1
#' # data(crm3)
#' # times<-crm3$times
#' # delta<-crm3$delta
#' # res <- CuMRes(times, delta, type.t = 2, 
#' #                   K = 100, length = .1, alpha = rep(1, 100  ), 
#' #                   beta = rep(1, 100),c.r = rep(50, 99), 
#' #                   iterations = 100, burn.in = 10, thinning = 1, type.c = 2)
#' 
#' 
#' @export CuMRes
CuMRes <-
  function(times, delta = rep(1, length(times)), type.t = 3, K = 5, utao = NULL,
           alpha = rep(0.01, K), beta = rep(0.01, K), 
           c.r = rep(1, (K - 1)),
           type.c = 4, epsilon = 1, c.nu = 1, a.eps = 0.1, b.eps = 0.1,
           a.mu = 0.01, b.mu = 0.01,
           iterations = 1000, burn.in = floor(iterations * 0.2), 
           thinning = 5, printtime = TRUE) {
    tInit <- proc.time()
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
    if (type.c %in% c( 2)) {
      if (length(c.r) != (K - 1)) {
        stop (c("Invalid argument: 'c.r' must have length, ", K - 1))
      }
      if (sum(abs(c.r - round(c.r)) > tol) != 0 || min(c.r) < 0) {
        stop ("Invalid argument: 'c.r' entries must be nonnegative integers.")
      }
    }
    if (type.c == 1 && sum(abs(c.r)) != (K-1) ) {
      c.r <- rep(0, K - 1)
      warning (c("'c.r' redefined as rep(0,", K - 1, ") because type.c = 1."))
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
    if (class(thinning) != "numeric") {
      stop ("Invalid argument: 'thinning' must be a logical value.")
    }
    if (thinning <= 0 || abs(thinning - round(thinning)) > tol 
        || thinning > 0.1 * iterations) {
      stop ("Invalid argument: 'thinning' must be a postitive integer smaller than 
            iterations * 0.10 = ", iterations * 0.1, ".")
    }
    if (printtime != TRUE && printtime != FALSE) {
      stop ("Invalid argument: 'printtime' must be a logical value.")
    }
    nm <- NM(times, delta, type.t, K, utao)
    n <- nm$n
    m <- nm$m
    tao <- nm$tao
    t.unc <- nm$t.unc
    acceptance.c <- 0
    if (type.c %in% c(3,4)) {
      c.r <- rep(5, (K - 1))
      Epsilon <- rep(NA, iterations)
    }
    cat(c("Iterating...", "\n"), sep = "")
    Lambda <- matrix(NA, nrow = iterations, ncol = K)
    U <- matrix(NA, nrow = iterations, ncol = K - 1)
    C <- matrix(NA, nrow = iterations, ncol = K - 1)
    lambda.r <- rep(0.1, K)
    Mu <- rep(NA, iterations)
    Z <- rep(NA, iterations) #iniciar vector de tiempo de quiebre
    Pi <- rep(NA, iterations) #iniciar vector de probabilidades
    k.star <- min(which(max(times[delta==1]) <= tao)) - 1 #k m?s grande donde hay al menos una observaci?n exacta
    z <- k.star # inicial para tiempo de quiebre
    pb <- dplyr::progress_estimated(iterations)
    for(j in seq_len(iterations)) {
      pb$tick()$print()
      u.r <- UpdU(alpha, beta, c.r, lambda.r)
      lambda.r <- CuUpdLambda(alpha, beta, c.r, u.r, n, m, z)
      mu <- rgamma(1, shape = a.mu + z - 1, rate = b.mu + 1) # simular mu
      aux.pi <- sum(lambda.r[seq_len(z)] * (tao[seq_len(z) + 1] - tao[seq_len(z)]))
      prop.pi <- exp(-aux.pi)
      z <- CuUpdZ(mu, m, lambda.r, k.star) # actualizar tiempo de quiebre
      if (type.c %in% c(3,4)) {
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
      Mu[j] <- mu
      Pi[j] <- prop.pi
      Z[j] <- z
      if (type.c %in% c(3,4)) Epsilon[j] <- epsilon
    }
    Lambda <- Lambda[seq(burn.in + 1, iterations, thinning), ]
    U <- U[seq(burn.in + 1, iterations, thinning), ]
    C <- C[seq(burn.in + 1, iterations, thinning), ]
    Mu <- Mu[seq(burn.in + 1, iterations, thinning)]
    Pi <- Pi[seq(burn.in + 1, iterations, thinning)]
    Z <- Z[seq(burn.in + 1, iterations, thinning)]
    Lambda <-  purrr::map_dfc(seq_len(ncol(Lambda)), ~ as.numeric((.x <= Z)*Lambda[,.x]))
    if (type.c %in% c(3,4)){ Epsilon <- Epsilon[seq(burn.in + 1, iterations, thinning)]}
    writeLines(c("","Done.", "\n", "Generating survival function estimates"))
    rows <- nrow(Lambda)
    s <- max(tao) * seq.int(0,100) / 100
    
    X <- as.matrix(unname(Lambda))
    pb <- dplyr::progress_estimated(length(s))
    S <- purrr::map(s, function(s = .x){
      pb$tick()$print()
      do.call(base::c, purrr::map(seq_len(rows),.f= ~exp(-sum((s > tao[-1]) * tao[-1] * X[.x,] + 
                                                         (s > tao[-length(tao)] & s <= tao[-1]) * s * X[.x,] -
                                                         (s > tao[-length(tao)]) * tao[-(length(tao))] * X[.x,])
      )))
    })
    
    cat(c("Done.", "\n"), sep = "")
    if (printtime) {
      cat(">>> Total processing time (sec.):\n")
      print(procTime <- proc.time() - tInit)
    }
    if(type.c %in% c(3,4)) {
      X = list(Lambda = Lambda, 
               U = U, C = C, Mu = Mu, Pi = Pi, Z = Z, Epsilon = Epsilon)} else { 
                 X = list(Lambda = Lambda, U = U, C = C, Mu = Mu, Pi = Pi, Z = Z)
               }
    X <- tibble::enframe(X) 
    out <- list(times = times, delta = delta, type.t = type.t, tao = tao, K = K, 
                t.unc = t.unc, iterations = rows, simulations = X, s = s,
                acceptance = acceptance.c/((K-1)*iterations),
                S = S)
    out <- tibble::enframe(out)
  }


