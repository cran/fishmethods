\name{agesurvcl}
\alias{agesurvcl}
\title{Age-Based Survival and Mortality Estimators for Cluster Sampling}
\description{
  Calculates the survival and mortality estimators of Jensen (1996) where net hauls are treated as samples
}
\usage{
agesurvcl(age = NULL, group = NULL, full = NULL, last = NULL)
}
\arguments{
  \item{age}{the vector of ages.  Each row represents the age of an individual.}
  \item{group}{the vector containing variable used to identify the sampling unit (e.g., haul). Identifier can be numeric or character.}
  \item{full}{the fully-recruited age.}
 \item{last}{the maximum age to include in the calculation. If not specified, the oldest age is used.}

}
\details{
  The individual age data are tabulated and subsetted based on \code{full} and \code{last}.  The calculations follow Jensen(1996).
 If only two ages are present, a warning message is generated.
}
\value{Matrix containing estimates of annual mortality (a), annual survival (S), and instantaneous total mortality (Z) and
associated standard errors.
}
\references{Jensen, A. L. 1996. \emph{Ratio estimation of mortality using catch curves}. Fisheries Research 27: 61-67.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{ \code{\link{agesurv}}}
\examples{
data(Jensen)
agesurvcl(age=Jensen$age,group=Jensen$group,full=0)
}
\keyword{misc}
