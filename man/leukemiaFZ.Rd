\name{leukemiaFZ}
\alias{leukemiaFZ}
\docType{data}
\title{
Survival Times of Patients with Leukemia (Feigl and Zelen)
}
\description{
Survival times of 33 patients with leukemia (Feigl and Zeelen, 1965). Times are measured in weeks from diagnosis. Reported covariates are white blood cell counts (WBC) and a binary variable AG that indicates a positive or negative test related to the white blood cell
characteristics. Three of the observations were censored. The data was taken from Lawless (2003).
}

\usage{data(leukemiaFZ)}

\format{
  A data frame with 33 observations on the following 4 variables.
  \describe{
    \item{\code{time}}{Weeks from diagnosis.}
    \item{\code{delta}}{Status indicator: 0=censored.}
    \item{\code{AG}}{Indicates a positive or negative test related to the white blood cell characteristics. (1=AG-positive, 2=AG-negative).}
    \item{\code{wbc}}{White blood cell counts in thousands (reported covariates).}
  }
}

\source{
Lawless, J.F. (2003). Statistical Models and Methods for Lifetime Data. Wiley: New Jersey.
}

\references{
Feigl, P. and Zelen, M. (1965). Estimation of Exponential Survival Probabilities with Concomitant Information. \emph{Biometrics} \strong{21}, 826-838.
}

\examples{
## Cox-Gamma Process Example 1
#  data(leukemiaFZ)
#  leukemia1 <- leukemiaFZ
#  leukemia1$wbc <- log(leukemiaFZ$wbc)
#  CGEX1 <- CGaMRes(data = leukemia1, K = 10, iterations = 100, thinning = 1)
}

\keyword{datasets}