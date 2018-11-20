\name{alkss}
\alias{alkss}
\title{Sample Size Determination for Age Subsampling}
\description{
	Calculates sample sizes for age subsampling assuming a two-stage random sampling design with
      proportional or fixed allocation.
}
\usage{
alkss(x, lss = NULL, cv = NULL, allocate = 1)
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
  \item{cv}{the desired coefficient of variation}
  \item{allocate}{the type of allocation: 1=proportional, 2=fixed.}
}
\details{
  If individual fish from catches are sampled randomly for lengths and then are further subsampled 
  for age structures, Quinn and Deriso (1999: pages 306-309) showed that sample sizes for age structures 
  can be determined for proportional (the number of fish aged is selected proportional to the length frequencies) and 
  fixed (a constant number are aged per length class) allocation assuming a two-stage random sampling design. 
  Sample sizes are determined based on the length frequency sample size, a specified 
  coefficient of variation, and proportional or fixed allocation.  The number of age classes is calculated internally.
See reference to Table 8.6, p. 312 in Quinn and Deriso.
  }
\value{
  \item{label }{list element containing the summary of input criteria}
  \item{n }{list element containing the sample size estimates for each age}
}
\references{Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{alkD}} \code{\link{alkprop}}}
\examples{
data(alkdata) 
alkss(alkdata,lss=1000,cv=0.25,allocate=1)
}
\keyword{misc}
