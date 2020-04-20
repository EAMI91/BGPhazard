\name{BMTKleinbook}
\alias{BMTKleinbook}
\docType{data}
\title{
Times to death or relapse (in days) for patients with bone marrow transplants for Hodgkin and non-Hodgkin lymphoma
}
\description{
Was collected on 43 bone marrow transplant patients at The Ohio State University Bone Marrow Transplant Unit. Details of this study can be found in Avalos et al. (1993).
}

\usage{data(BMTKleinbook)}

\format{
  A data frame with 43 observations containing:
  \describe{
    \item{\code{times}}{time to death or relapse in days}
    \item{\code{delta}}{Status indicator: 1 = death or relapse; 0 = otherwise}
    \item{\code{tTransplant}}{Allogeneic transplant from an HLA match sibling donor (1) or an autogeneic transplant (0)}
    \item{\code{hodgkin}}{Hodgkin disease (1), or non-Hodgkin lymphoma (0)}
    \item{\code{karnofsky}}{The pretransplant Karnofsky score}
    \item{\code{waiting}}{Waiting time to transplant}
  }
}

\source{
Klein, J. P., and Moeschberger, M. L. (2003). Survival analysis: techniques for censored and truncated data. Springer Science & Business Media. 
}

\references{
Copelan, E. A., Biggs, J. C., Thompson, J. M., Crilley, P., Szer, J., Klein, J. P., Kapoor, N.,
Avalos, B. R., Cunningham, I., Atkinson, K., Downs, K., Harmon, G. S., Daly, M. B., Brodsky, I., Bulova, S. I., and Tutschka, P. J. Treatment for Acute Myelocytic Leukemia with Allogeneic Bone Marrow Transplantation Following Preparation with Bu/Cy. Blood 78 (1991): 838-843.
}

\examples{
## Cox Cure Gama Process Example 1
# data(BMTKleinbook)
    # res <- CCuMRes(BMTKleinbook, covs.x = c("tTransplant","hodgkin","karnofsky","waiting"),
    #                covs.y = c("tTransplant","hodgkin","karnofsky","waiting"),
    #                        type.t = 2, K = 72, length = 30,
    #                        alpha = rep(2,72), beta = rep(2,72), c.r = rep(50, 71), type.c = 2,
    #                        var.delta.str = .1, var.theta.str = 1,
    #                        var.delta.ini = 100, var.theta.ini = 100,
    #                        iterations = 100, burn.in = 10, thinning = 1)
}

\keyword{datasets}