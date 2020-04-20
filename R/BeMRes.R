#' Markov Beta Model
#' 
#' Posterior inference for the Bayesian non-parametric Markov beta model for discrete
#' survival times.
#' 
#' Computes the Gibbs sampler given by the full conditional distributions of u
#' and Pi (Nieto-Barajas & Walker, 2002) and arranges the resulting Markov
#' chain into a tibble which can be used to obtain posterior summaries.
#' 
#' @param times Numeric positive vector. Failure times.
#' @param delta Logical vector. Status indicator. \code{TRUE} (1) indicates
#' exact lifetime is known, \code{FALSE} (0) indicates that the corresponding
#' failure time is right censored.
#' @param alpha Nonnegative vector. Small entries are recommended in order to
#' specify a non-informative prior distribution.
#' @param beta Nonnegative vector. Small entries are recommended in order to
#' specify a non-informative prior distribution.
#' @param c.r Nonnegative vector. The higher the entries, the higher the
#' correlation of two consecutive failure times.
#' @param a.eps Numeric. Shape parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param b.eps Numeric. Scale parameter for the prior gamma distribution of
#' epsilon when \code{type.c = 4}.
#' @param type.c Integer. 1=defines \code{c.r} as a zero-entry vector; 2=lets
#' the user define \code{c.r} freely; 3=assigns \code{c.r} an
#' exponential prior distribution with mean \code{epsilon}; 4=assigns \code{c.r} 
#' an exponential hierarchical distribution with mean \code{epsilon} which in turn has a
#' a Ga(a.eps, b.eps) distribution.
#' @param epsilon Double. Mean of the exponential distribution assigned to
#' \code{c.r}
#' @param iterations Integer. Number of iterations including the \code{burn.in}
#' and \code{thining} to be computed for the Markov chain.
#' @param burn.in Integer. Length of the burn-in period for the Markov chain.
#' @param thinning Integer. Factor by which the chain will be thinned. Thinning
#' the Markov chain is to reduces autocorrelation.
#' @param printtime Logical. If \code{TRUE}, prints out the execution time.
#' @note It is recommended to verify chain's stationarity. This can be done by
#' checking each partition element individually. See \link{BePlotDiag}.
#' @seealso \link{BePlotDiag}, \link{BePloth}
#' @references - Nieto-Barajas, L. E. & Walker, S. G. (2002). Markov beta and
#' gamma processes for modelling hazard rates. \emph{Scandinavian Journal of
#' Statistics} \strong{29}: 413-424.
#' @examples
#' 
#' 
#' 
#' ## Simulations may be time intensive. Be patient.
#' 
#' ## Example 1
#' #  data(psych)
#' #  timesP <- psych$time
#' #  deltaP <- psych$death
#' #  BEX1 <- BeMRes(timesP, deltaP, iterations = 3000, burn.in = 300, thinning = 1)
#' 
#' ## Example 2
#' #  data(gehan)
#' #  timesG <- gehan$time[gehan$treat == "control"]
#' #  deltaG <- gehan$cens[gehan$treat == "control"]
#' #  BEX2 <- BeMRes(timesG, deltaG, type.c = 2, c.r = rep(50, 22))
#' 
#' 
#' 
#' @export BeMRes
BeMRes <-
  function(times, delta = rep(1, length(times)), alpha = rep(0.0001, K), 
           beta = rep(0.0001, K), c.r = rep(0, K-1), 
           a.eps = 0.1, b.eps = 0.1, type.c = 4, 
           epsilon =  1, iterations = 2000, 
           burn.in = floor(iterations * 0.2), thinning = 5, printtime = TRUE) {  
    tInit <- proc.time()
    K <- max(times)
    tol = .Machine$double.eps ^ 0.5
    if (min(times) < 0) {
      stop ("Invalid argument: 'times' must be a nonnegative integer vector.")
    }
    if (max(abs(times - round(times)) > tol) == 1) {
      stop ("Invalid argument: 'times' must be a nonnegative integer vector.")
    } 
    if (min((delta  ==  0) + (delta  ==  1 )) == 0) {
      stop ("Invalid argument: 'delta' must have 0 - 1 entries.")
    }
    if (length(times) != length(delta)) {
      stop ("Invalid argument: 'times' and 'delta' must have same length.")
    }
    if (length(alpha) != K || length(beta) != K) {
      stop (c("Invalid argument: 'alpha', 'beta', must have length "), K, ("."))
    }
    if (min(c(alpha, beta)) < 0) {
      stop ("Invalid argument: 'alpha' and 'beta' must have nonnegative entries.")
    } 
    if (abs(type.c - round(type.c)) > tol || type.c < 1 || type.c > 4) {
      stop ("Invalid argument: 'type.c' must be an integer between 1 and 4.")
    }
    if (type.c == 2) {
      if (length(c.r) != (K - 1)) {
        stop (c("Invalid argument: 'c.r' must have length, ", K - 1, "."))
      }
      if (max(abs(c.r - round(c.r)) > tol) == 1 || min(c.r) < 0) {
        stop ("Invalid argument: 'c.r' entries must be positive integers.")
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
        || burn.in > iterations * 0.9) {
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
    nm <- BeNM(times, delta)
    n <- nm$n
    m <- nm$m
    tao <- nm$tao
    t.unc <- nm$t.unc
    if (type.c == 3 || type.c == 4) {
      c.r <- rep(5, (K - 1))
    }
    if (type.c == 4) {
      Epsilon <- rep(NA, iterations)
    }
    cat(c("Iterating...", "\n"), sep = "")
    PI <- matrix(NA, nrow = iterations, ncol = K)
    U <- matrix(NA, nrow = iterations, ncol = (K - 1))
    C <- matrix(NA, nrow = iterations, ncol = (K - 1))
    Pi.r <- rep(0.1, K)
    pb <- dplyr::progress_estimated(iterations)
    for(j in seq_len(iterations)) {
      pb$tick()$print()
      u.r <- BeUpdU(alpha, beta, c.r, Pi.r)
      Pi.r <- UpdPi(alpha, beta, c.r, u.r, n, m)
      if (type.c == 3 || type.c == 4) {
        if (type.c == 4) {
          epsilon <- rgamma(1, shape = a.eps + K, scale = 1 / (b.eps + sum(c.r)))
        }
        c.r <- BeUpdC(alpha, beta, Pi.r, u.r, epsilon)
      }
      PI[j, ] <- Pi.r
      U[j, ] <- u.r
      C[j, ] <- c.r
      if(type.c == 4) Epsilon[j] <- epsilon
    }
    PI <- PI[seq(burn.in + 1, iterations, thinning), ]
    U <- U[seq(burn.in + 1, iterations, thinning), ]
    C <- C[seq(burn.in + 1, iterations, thinning), ]
    if (type.c == 4) Epsilon <- Epsilon[seq(burn.in + 1, iterations, thinning)]
    cat(c("Done.", "\n", "Generating survival function estimates.", "\n"),
        sep = "")
    rows <- nrow(PI)
    s <- seq_len(K)
    S <- tibble::as_tibble(do.call(rbind,unname(purrr::map_dfr(purrr::map(tibble::as_tibble(t(PI)) ,~purrr::accumulate(.x,`+`)),~exp(-.x)))))
    if (printtime) {
      cat(">>> Total processing time (sec.):\n")
      print(procTime <- proc.time() - tInit)
    }
    if(type.c == 4) {
      X = tibble::enframe(list(PI = PI, 
               U = U, C = C, Epsilon = Epsilon))
    } else { 
      X = tibble::enframe(list(PI = PI, U = U, C = C))
    }
    out <- tibble::enframe(list(times = times, delta = delta, tao = tao, K = K, t.unc = t.unc,
                iterations = rows, simulations = X, s = s, S = S))
    out
  }
