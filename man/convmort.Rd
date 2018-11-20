\name{convmort}
\alias{convmort}
\title{Conversion of Mortality Rates}
\description{
  Convert instantaneous fishing mortality rate (F) to annual exploitation rate (mu) and vice versa for
Type I and II fisheries.
}
\usage{
convmort(value = NULL, fromto = 1, type = 2, M = NULL)
}
\arguments{
  \item{value}{mortality rate}
  \item{fromto}{conversion direction: 1=from F to mu; 2 = from mu to F.  Default is 1.}
  \item{type}{type of fishery following Ricker (1975): 1=Type I; 2=Type II. Default is 2.}
  \item{M}{natural mortality rate (for Type II fishery)}
}
\details{
  Equations 1.6 and 1.11 of Ricker (1975) are used.
}
\value{A vector of the same length as \code{value} containing the converted values.
}
\references{Ricker, W. E. 1975. Computation and interpretation of biological statistics of fish populations.
Bull. Fish. Res. Board. Can. 191: 382 p.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
convmort(0.3,fromto=1,type=2,M=0.15)
}
\keyword{misc}

