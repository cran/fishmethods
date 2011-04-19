\name{alkprop}
\alias{alkprop}
\title{Age-Length Key Proportions-At-Age}
\description{
	Calculates proportions-at-age and standard errors from an age-length key assuming a two-stage random sampling design.
}
\usage{
alkprop(x)
}
\arguments{
  \item{x}{a data frame containing an age-length key (similar to Table 8.3 on page 307 of Quinn and Deriso (1999)).  
   The first column must contain the length intervals as single numeric labels (no ranges), the second column must 
   contain the number of samples within each length interval (Ll in Q & D), and the third and remaining columns must 
   contain the number of samples for each age class within each length interval (one column for each age class).
   Column labels are not necessary but are helpful. Columns l and Al in Table 8.3 should not be 
   included. Empty cells must contain zeros.
  }
}
\details{
  If individual fish from catches are sampled randomly for lengths and then are further subsampled 
  for age structures, Quinn and Deriso (1999: pages 304-305) showed that the proportions of fish in each age class
  and corresponding standard errors can be calculated assuming a two-stage random sampling design. See reference to Table 8.4,
page 308 in Quinn and Deriso.
  }
\value{
  \item{results}{list element containing a table of proportions, standard errors, and coefficients of variation
  for each age class.}}
\references{Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{alkD}} \code{\link{alkss}}}
\examples{
data(alkdata) 
alkprop(alkdata)
}
\keyword{misc}
