\name{cpuekapp}
\alias{cpuekapp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trawl survey based abundance estimation using data sets with unusually large catches
}
\description{
Calculates the mean cpue after replacing unusually large catches with
expected values using the method of Kappenman (1999) 
}
\usage{
cpuekapp(x = NULL, nlarge = NULL, absdif = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of non-zero trawl catch data.
}
  \item{nlarge}{the number of values considered unusually large.
}
  \item{absdif}{convergence tolerance}
}
\details{
Use function \code{gap} to choose the number of unusually large values.

}
\value{
 \item{kappmean}{list element containing new arithmetic mean.}
 \item{expectations }{list element containing the original observation(s)
  and expected order statistic(s).}
}
\references{
Kappenman, R. F. 1999. Trawl survey based abundance estimation using data sets with 
unusually large catches. ICES Journal of Marine Science. 56: 28-35.
}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}
}
\seealso{\code{\link{gap}}}
\examples{
 \dontrun{
  ## Data from Table 1 in Kappenman (1999)
  data(kappenman)
  cpuekapp(kappenman$cpue,1)
 }
}
\keyword{misc}
