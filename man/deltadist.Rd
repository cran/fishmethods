\name{deltadist}
\alias{deltadist}
\title{Delta Distribution Mean and Variance Estimators}
\description{
  Calculates the mean and variance of a catch series based on the delta distribution described in 
Pennington (1983).
}
\usage{
deltadist(x = NULL)
}
\arguments{
  \item{x}{vector of catch values, one record for each haul.  Include zero and nonzero catches. Missing values are deleted prior to estimation.}
}
\details{
Data from marine resources surveys usually contain a large proportion of hauls
with no catches.  Use of the delta-distribution can lead to more efficient estimators of 
the mean and variance because zeros are treated separately.  The methods used here to calculate
the delta distribution mean and variance are given in Pennington (1983).  
}
\value{
  vector containing the delta mean and associated variance.
  
}
\references{Pennington, M. 1983. Efficient estimators of abundance for fish and plankton
surveys. Biometrics 39: 281-286.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
data(catch)
deltadist(catch$value) 
}            
\keyword{misc}