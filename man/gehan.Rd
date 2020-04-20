\name{gehan}
\alias{gehan}
\docType{data}
\title{
Times of Remission of Leukemia Patients (Gehan)
}
\description{
Freireich et al. (1963) report the results of a clinical trial of a drug 6-mercaptopurine (6-MP) versus a placebo (control) in 42 children with acute leukemia. The trial was conducted at 11 American hospitals. The trial was conducted by matching pairs of patients at a given hospital by remission status (complete or partial) and randomizing within the pair to either a 6-MP or placebo maintenance therapy. Patients were followed until their leukemia returned (relapse) or until the end of the study (in weeks). The data was taken from Klein & Moeschberger (2003) and is contained in the \code{MASS} package.
}

\usage{data(gehan)}

\format{
  A data frame with 42 observations containing:
  \describe{
    \item{\code{pair}}{Pair index.}
    \item{\code{time}}{Remission time (weeks).}
    \item{\code{cens}}{Status: 0=censored.}
    \item{\code{treat}}{Treatment: control or 6-MP.}
  }
}

\source{
Klein, J. P., & Moeschberger, M. L. (2003). Survival analysis: techniques for censored and truncated data. Springer Science & Business Media. 
}

\references{
Freireich, E. J., et al. (1963). The effect of 6-mercaptopurine on the duration of steroid-induced remissions in acute leukemia: A model for evaluation of other potentially useful therapy. \emph{Blood}, \strong{21(6)}, 699-716.
}

\examples{
## Gamma Process Example 1 
#  data(gehan)
#  timesG <- gehan$time[gehan$treat == "6-MP"]
#  deltaG <- gehan$cens[gehan$treat == "6-MP"]
#  GEX1 <- GaMRes(timesG, deltaG, K = 8, iterations = 3000)
#  GaPloth(GEX1)
}

\keyword{datasets}
