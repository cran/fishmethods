\name{alk}
\alias{alk}
\title{Create An Age-Length Key}
\description{
	Creates an age-length key in numbers or proportions-at-age per size.
}
\usage{
alk(age=NULL,size=NULL,binsize=NULL,type=1)
}
\arguments{
  \item{age}{a vector of individual age data.}
  \item{size}{a vector of individual size data.}
  \item{binsize}{size of the length class (e.g., 5-cm, 10, cm, etc.) used to group size data.
       The formula used to create bins is \eqn{trunc(size/binsize)*binsize+binsize/2}. 
       If use of the raw length classes is desired, then \code{binsize=0}.}
  \item{type}{If \code{type}=1, numbers-at-age per size are produced.  This format is used in
    functions \code{alkprop}, \code{alkss}, and \code{alkD}. If \code{type}=2,
    proportions-at-age per size are produced.}
}
\details{
 Create age-length keys with either numbers-at-age per size class. Records with missing size values are deleted prior to calculation.
     Missing ages are allowed.
  }
\value{A table of size, total numbers at size, and numbers (or proportions)-at-age per size class.}
\references{Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{alkD}} \code{\link{alkss}} \code{\link{alkprop}}}
\examples{
 \dontrun{
  data(pinfish) 
  with(pinfish,alk(age=round(age,0),size=sl,binsize=10))
 }
}
\keyword{misc}
