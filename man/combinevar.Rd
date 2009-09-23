\name{combinevar}
\alias{combinevar}
\title{Combining Mean and Variances from Multiple Samples}
\description{
 This function takes multiple mean and sample variance estimates and combines them. 
}
\usage{
combinevar(xbar = NULL, s_squared = NULL, n = NULL)
}
\arguments{
  \item{xbar}{vector of means}
  \item{s_squared}{vector of sample variances}
  \item{n}{vector of number of observations}
}
\details{
If a Monte Carlo simulation is run over 1000 loops and then again over another 1000 loops, one may wish to update 
the mean and variance from the first 1000 loops with the second set of simulation results.
}
\value{Vector containing the combined mean and sample variance.
}
\author{John M. Hoenig, Virginia Institute of Marine Science \email{hoenig@vims.edu}}
\examples{
xbar <- c(5,5)
s<-c(2,4)
n <- c(10,10)
combinevar(xbar,s,n)
}
\keyword{misc}

