\name{mort.al}
\alias{mort.al}
\title{
 Estimation of Mortality using Times-At-Large Data from Tagging
}
\description{
 Calculates total instantaneous (Z), natural mortality (M) and/or fishing mortality (F) using times-at-large data
and methods of Gulland (1955) and McGarvey et al. (2009).
}
\usage{
mort.al(relyr = NULL, tal = NULL, N = NULL, method = c(1, 2, 3), 
np = 0, stper = NULL, nboot = 500)
}
\arguments{
  \item{relyr}{a vector of release year (or cohort) for individual times-at-large observations.}
  \item{tal}{a vector of individual times-at-large observations.}
  \item{N}{a vector of number of releases for each release year (or cohort). Each individual observation
	from a release year should have the same N value.}
  \item{method}{1 = McGarvey et al., 2 = Gulland.  Default is all (i.e., c(1,2)).}
  \item{np}{the number of periods over which to combine data to make period estimates of mortality. Set 
   np=0 to estimate mortality for each release year.}
  \item{stper}{vector of year values representing the beginning year of each period over which 
	to estimate mortality. The first year in c() must always be the first release year.}
  \item{nboot}{the number of resamples for the Gulland method.}
}
\details{
  The methods of Gulland (1955) and McGarvey et al (2009) are used to estimate Z, F and M (depending
on the method) from tagging times-at-large data.  For the Gulland method,  the standard error of the Z, M, and F 
estimates are made using a parametric bootstrap method similar to Tanaka (2006). When
 periods are specified, period-specific mortality estimates and standard errors are derived by
averaging release-year-specific mortality estimates. The standard errors are calculated by taking the square-root of 
the averaged variances of the estimates. To combine data over all years prior to estimation, change all relyr within 
a period to the same year value.

}
\value{dataframe containing the M, F and Z estimates and associated standard errors by period.}
\references{
Gulland, J. A. 1955. On the estimation of population parameters from marked members. Biometrika 42: 269-270.

McGarvey, R., J. M. Matthews, and J. E. Feenstra. 2009.  Estimating mortality from times-at-large: testing accuracy and precision 
using simulated single tag-recovery data. ICES Journal of Marine Science 66: 573-581.

Tanaka, E. 2006. Simultaneous estimation of instantaneous mortality coefficients and rate of effective survivors to number of released fish 
using multiple sets of tagging experiments. Fisheries Science 72: 710-718.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
 \dontrun{
  data(tanaka)
  mort.al(relyr = tanaka$relyr, tal = tanaka$tal, N = tanaka$N)
 }
}
\keyword{misc}