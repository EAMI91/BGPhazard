% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GaMRes.R
\name{GaMRes}
\alias{GaMRes}
\title{Markov Gamma Model}
\usage{
GaMRes(
  times,
  delta = rep(1, length(times)),
  type.t = 3,
  K = 5,
  utao = NULL,
  alpha = rep(0.01, K),
  beta = rep(0.01, K),
  c.r = rep(1, (K - 1)),
  c.nu = 1,
  a.eps = 0.1,
  b.eps = 0.1,
  type.c = 4,
  epsilon = 1,
  iterations = 1000,
  burn.in = floor(iterations * 0.2),
  thinning = 5,
  printtime = TRUE
)
}
\arguments{
\item{times}{Numeric positive vector. Failure times.}

\item{delta}{Logical vector. Status indicator. \code{TRUE} (1) indicates
exact lifetime is known, \code{FALSE} (0) indicates that the corresponding
failure time is right censored.}

\item{type.t}{Integer. 1=computes uniformly-dense intervals; 2=
partition arbitrarily defined by the user with parameter utao and 3=same length intervals.}

\item{K}{Integer. Partition length for the hazard function if
\code{type.t}=1 or \code{type.t}=3.}

\item{utao}{vector. Partition specified by the user when type.t = 2. The first value of 
the vector has to be 0 and the last one the maximum observed time, either censored or uncensored.}

\item{alpha}{Nonnegative entry vector. Small entries are recommended in
order to specify a non-informative prior distribution.}

\item{beta}{Nonnegative entry vector. Small entries are recommended in order
to specify a non-informative prior distribution.}

\item{c.r}{Nonnegative vector. The higher the entries, the higher the correlation 
of two consecutive intervals.}

\item{c.nu}{Tuning parameter for the proposal distribution for c.}

\item{a.eps}{Numeric. Shape parameter for the prior gamma distribution of
epsilon when \code{type.c = 4}.}

\item{b.eps}{Numeric. Scale parameter for the prior gamma distribution of
epsilon when \code{type.c = 4}.}

\item{type.c}{1=assigns \code{c.r} a zero-entry vector; 2=lets the user
define \code{c.r} freely; 3=assigns \code{c.r} an exponential prior
distribution with mean 1; 4=assigns \code{c.r} an exponential hierarchical
distribution with mean \code{epsilon} which in turn has a Ga(a.eps, b.eps)
distribution.}

\item{epsilon}{Double. Mean of the exponential distribution assigned to
\code{c.r} when \code{type.c = 3}}

\item{iterations}{Integer. Number of iterations including the \code{burn.in}
to be computed for the Markov chain.}

\item{burn.in}{Integer. Length of the burn-in period for the Markov chain.}

\item{thinning}{Integer. Factor by which the chain will be thinned. Thinning
the Markov chain is to reducec autocorrelation.}

\item{printtime}{Logical. If \code{TRUE}, prints out the execution time.}
}
\description{
Computes the Gibbs sampler given by the full conditional distributions of U,
Lambda, C and Epsilon (Nieto-Barajas & Walker, 2002) and arranges the
resulting Markov chain into a tibble which can be used to obtain posterior
summaries.
}
\details{
Posterior inference for the Bayesian non-parametric Markov gamma model in
survival analysis.
}
\examples{


## Simulations may be time intensive. Be patient.

## Example 1 
data(gehan) 
timesG <- gehan$time[gehan$treat == "6-MP"] 
deltaG <- gehan$cens[gehan$treat == "6-MP"] 
 GEX1 <- GaMRes(timesG, deltaG, K = 8, iterations = 3000)

## Example 2 
data(leukemiaFZ) 
timesFZ <- leukemiaFZ$time 
deltaFZ <- leukemiaFZ$delta 
GEX2 <- GaMRes(timesFZ, deltaFZ, type.c = 4)



}
