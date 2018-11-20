\name{alkD}
\alias{alkD}
\title{Sample Size Determination for Age Subsampling Using the D statistic}
\description{
Calculates the D statistic (sqrt of accumulated variance among ages; Lai 1987) for a range of age sample sizes using data 
from an age-length key.  Assumes a two-stage random sampling design with proportional or fixed allocation.
}
\usage{
alkD(x, lss = NULL, minss = NULL, maxss = NULL, sampint = NULL, 
    allocate = 1)
}
\arguments{
  \item{x}{a data frame containing an age-length key (similar to Table 8.3 on page 307 of Quinn and Deriso (1999)).  
   The first column must contain the length intervals as numeric labels (no ranges), the second column must 
   contain the number of samples within each length interval (Ll in Q & D), and the third and remaining columns must 
   contain the number of samples for each age class within each length interval (one column for each age class).
   Column labels are not necessary but are helpful. Columns l and Al in Table 8.3 should not be 
   included. Empty cells must contain zeros.
}
  \item{lss}{the sample size for length frequency}
  \item{minss}{the minimum age sample size}
  \item{maxss}{the maximum age sample size.  Value can not be larger than the sample size for the length frequency(\code{lss})}
  \item{sampint}{the sample size interval}
  \item{allocate}{the type of allocation: 1=proportional, 2=fixed.}
}
\details{
  Following Quinn and Deriso (1999:pages 308-309), the function calculates the D statistic (sqrt of 
accumulated variance among ages; Lai 1987) for a range of age sample sizes defined by \code{minss}, \code{maxss}, and 
\code{sampint} at a
given length sample size \code{lss}.  The size of an age sample at a desired level of D can be obtained by the comparison. See reference
to Table 8.8, p. 314 in Quinn and Deriso.
}
\value{
  \item{label }{list element containing the summary of input criteria}
  \item{comp2 }{list element containing the D statistic for each age sample size given lss}
 }
\references{Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages

Lai, H.L. 1987. Optimum allocation for estimating age composition using age-length keys. U.S. Fish. Bull. 85:179-185}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{alkss}} \code{\link{alkprop}}}
\examples{
data(alkdata)
alkD(alkdata,lss=1000,minss=25,maxss=1000,sampint=20,allocate=1)
}
\keyword{misc}
