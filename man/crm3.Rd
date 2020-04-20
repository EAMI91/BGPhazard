\name{crm3}
\alias{crm3}
\docType{data}
\title{
Simulated data from the proportional hazards cure rate model of Yakolev & Tsodikov (1996) by taking a baseline density with a bounded support.
}
\description{
Triangular distribution Tri(a, c, b) as the baseline density, which puts a probability of one to the interval [a, b] and the mode at c. 
}
\usage{data("crm3")}
\format{
  A data frame with 100 observations with the following 2 variables.
  \describe{
    \item{\code{times}}{Simulated time}
    \item{\code{delta}}{Simulated censoring}
  }
}
\details{
In particular we took, a = 0, c = 1 and b = 4. The censoring time was independently generated from a uniform distribution to yield a 30\% censoring rate. Sample size n = 100 and the cure proportion exp\{-theta\}=0.2.
}
\source{
Nieto-Barajas, L. E., & Yin, G. (2008)
}

\examples{
## Cure Gama model Example 1
# data(crm3)
# times<-crm3$times
# delta<-crm3$delta
# res <- CuMRes(times, delta, type.t = 2, 
#                   K = 100, length = .1, alpha = rep(1, 100  ), 
#                   beta = rep(1, 100),c.r = rep(50, 99), 
#                   iterations = 100, burn.in = 10, thinning = 1, type.c = 2)
}
\keyword{datasets}
