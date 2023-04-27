\name{simulus}
\alias{simulus}
\docType{data}
\title{
Age and size data for the growth_sel function
}
\description{
 Age and size data were derived via simulation.
}
\usage{data(simulus)}
\format{
  A data frame with 1000 observations on the following 6 variables.
  \describe{
    \item{\code{age}}{a numeric vector of ages}
    \item{\code{size}}{a numeric vector of body size}
    \item{\code{weights}}{a numeric vector of observation weights for the likelihood function.}
    \item{\code{minlimit}}{a numeric vector of the minimum size limit.}
    \item{\code{maxlimit}}{a numeric vector of the maximum size limit.}
    \item{\code{minmax}}{a numeric vector indicating to which likelihood component (1=minimum, 2=maximum) each row observation is assigned. }
  }
}

\source{
Amy M. Schueller, National Marine Fisheries Service, Beaufort, NC \email{amy.schueller@noaa.gov}
}
\keyword{datasets}
