\name{psych}
\alias{psych}
\docType{data}
\title{
Death Times of Psychiatric Patients (Woolson)
}
\description{
Woolson (1981) has reported survival data on 26 psychiatric inpatients admitted to the University of Iowa hospitals during the years 1935-1948. This sample is part of a larger study of psychiatric inpatients discussed by Tsuang and Woolson (1977) and it is contained in the \code{KMsurv} package.
}

\usage{data(psych)}

\format{
  A data frame with 26 observations containing:
  \describe{
    \item{\code{sex}}{Patient sex: 1=male, 2=female.}
    \item{\code{age}}{Age at first admission to the hospital.}
    \item{\code{time}}{Number of years of follow-up.}
    \item{\code{death}}{Patient status at the follow-up time: 0=alive, 1=dead.}
  }
}

\source{
Klein, J. P., and Moeschberger, M. L. (2003). Survival analysis: techniques for censored and truncated data. Springer Science & Business Media. 
}

\references{
Tsuang, M. T. and Woolson, R. F. (1977). Mortality in Patients with Schizophrenia, Mania and
Depression. \emph{British Journal of Psychiatry}, \strong{130}: 162-166.

Woolson, R. F. (1981). Rank Tests and a One-Sample Log Rank Test for Comparing Observed
Survival Data to a Standard Population. \emph{Biometrics} \strong{37}: 687-696.
}

\examples{
## Beta Process Example 1
## Example 1
#  data(psych)
#  timesP <- psych$time
#  deltaP <- psych$death
#  BEX1 <- BeMRes(timesP, deltaP, iterations = 3000, burn.in = 300, thinning = 1)
}

\keyword{datasets}