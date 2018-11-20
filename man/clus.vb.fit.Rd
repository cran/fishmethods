\name{clus.vb.fit}
\alias{clus.vb.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Von Bertalanffy growth equation to clustered data via bootstrapping }
\description{
Fits the von Bertalanffy growth equation to clustered length and age  
by using nonlinear least-squares and by bootstrapping clusters}
\usage{clus.vb.fit(len=NULL,age=NULL,cluster=NULL,nboot=1000,
                    control=list(maxiter=10000,
                    minFactor=1/1024,tol=1e-5))
}
\arguments{
  \item{len}{vector of lengths of individual fish}
  \item{age}{vector of ages of individual fish}
  \item{cluster}{haul or cluster membership identifier}
  \item{nboot}{number of bootstrap samples}
  \item{control}{see \code{control} under function \emph{nls}.}
}
\details{
A standard von Bertalanffy growth curve is fitted to length-at-age data of each  \emph{nboot} sample of clusters
 using nonlinear least-squares (function \emph{nls}).  The estimated parameters and their associated standard deviations
  are calculated as the mean and variance of the estimates from all bootstrap fits.  The average parameter correlation matrix
from all bootstrap fits is provided as well.
}
\value{
 List containing a summary of successful model fits and parameter estimates, standard deviation and 
  95 percent confidence intervals, and the average correlation matrix.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\examples{
  \dontrun{
	data(pinfish)
	with(pinfish,clus.vb.fit(len=sl,age=age,cluster=field_no,nboot=100))
  }
}
\keyword{misc}
