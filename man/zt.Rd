\name{zt}
\alias{zt}
\title{Z-transform or center a time series}
\description{
  Z-transforms observations of a time series or centers observations of a time series to the mean.
}
\usage{
zt(x = NULL, ctype = 1)
}
\arguments{
  \item{x}{vector of observations. Missing values are allowed.}
  \item{ctype}{the type of transformation.  1 = Z transform ((x - mean x)/ sd x); 2 = center (x - mean x). Default = 1}
}
\details{
 Z-transforms observations of a time series or centers observations of a time series to the mean.
}
\value{
  vector containing the transformed time series.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
data(wolffish)
zt(wolffish$index) 
}            
\keyword{misc}