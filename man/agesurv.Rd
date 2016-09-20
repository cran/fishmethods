\name{agesurv}
\alias{agesurv}
\title{Age-based Survival Estimators}
\description{
  Calculates annual survival (S) and instantaneous total mortality rates (Z) from age frequency by using
catch curve (ln(numbers) versus age in linear regression), Heincke, Chapman-Robson and GLMER methods.
}
\usage{
agesurv(age = NULL, full = NULL, last = NULL, 
estimate = c("s", "z"), method = c("cc", "he", "cr", "crcb","mil","mb"))
}
\arguments{
  \item{age}{ the vector of ages.  Each row represents the age of an individual.}
  \item{full}{the fully-recruited age}
  \item{last}{the maximum age to include in the calculation. If not specified, the oldest age is used.}
  \item{estimate}{argument to select estimate type: "s" for annual survival, "z" for instantaneous total mortality.
	Default is both.}
  \item{method}{argument to select the estimation method: "cc" for catch curve, "he" for Heincke,
       "cr" for Chapman-Robson, "crcb" for Chapman-Robson Z estimate bias-correction (Seber p. 418) and 
       SEZ corrected for over-dispersion (Smith et al., 2012), "mil" for Millar (2015) random-intercept
	poisson log-linear mixed model estimator, "mb" for Maceine-Bettoli weighted regression.  Default is all.}
}
\details{The individual age data are tabulated and subsetted based on the \code{full} and \code{last} arguments.  
  Most calculations follow descriptions in Seber(1982), pages 414-418. If only two ages are present, a warning message 
  is generated and the catch curve method is not calculated. Plus groups are not allowed.
}
\value{
  \item{results}{list element containing table of parameters and standard errors.}
  \item{data}{list element containing the age frequency data used in the analysis.}
}
\references{
Seber, G. A. F. 1982. The Estimation of Animal Abundance and Related Parameters, Second Edition. The Blackburn Press,
Caldwell, New Jersey. 654 pages. 

Maceina, M. J. and P. W. Bettoli. 1998. Variation in largemouth bass recruitment in four mainstream impoundments of the Tennessee River. N. Am. J. Fish. Manage. 18: 990-1003.

Millar, R. B. 2015. A better estimator of mortality rate from age-frequency data. Can. J. Fish. Aquat.
Sci. 72: 364-375.

Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages.

Smith, M. W. and 5 others. 2012. Recommendations for catch-curve analysis. N. Am. J. Fish. Manage. 32: 956-967.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{data(rockbass)
agesurv(age=rockbass$age,full=6)}
\keyword{misc}
