\name{agesurv}
\alias{agesurv}
\title{Age-based Survival Estimators}
\description{
  Calculates annual survival (S) and instantaneous total mortality rates (Z) from age frequency by using
linear regression (standard and weighted), Heincke, Chapman-Robson, Poisson GLM and GLMER methods.
}
\usage{
agesurv(type=1, age=NULL, number=NULL, full=NULL, last=NULL, estimate=c("s","z"),
method=c("lr","he","cr","crcb","ripois","wlr","pois"), sign.est=3, sign.se=3, 
 glmer.control=glmerControl(optCtrl=list(maxfun=10000),optimizer="bobyqa"))
}

\arguments{
  \item{type}{ the format of data. 1 = a single vector, each row represents the age of an individual (default), 
        2 = summarized, one column of age and one column of numbers-at-age.}

  \item{age}{ the vector of ages.}
  \item{number}{ if \code{type} = 2, a vector of numbers-at-age matching the length of the age vector.}
  \item{full}{the fully-recruited age}
  \item{last}{the maximum age to include in the calculation. If not specified, the oldest age is used.}
  \item{estimate}{argument to select estimate type: "s" for annual survival, "z" for instantaneous total mortality.
	Default is both.}
  \item{method}{argument to select the estimation method: "lr" for standard linear regression, "he" for Heincke,
       "cr" for Chapman-Robson, "crcb" for Chapman-Robson Z estimate with bias-correction (Seber p. 418) and 
        over-dispersion correction (Smith et al., 2012), "ripois" for Millar (2015) random-intercept
	Poisson mixed model estimator, "wlr" for Maceine-Bettoli weighted regression,
        "pois" for Poisson generalized linear model with overdispersion correction.
         Default is all.}
  \item{sign.est}{significant digits for survival estimates.}
  \item{sign.se}{significant digits for standard error of survival estimates.}

  \item{glmer.control}{controls for function \code{glmer} used in the random-intercept Poisson mixed model. See
       \code{glmerControl}.}
}
\details{If \code{type} = 1, the individual age data are tabulated. The age data are then subsetted based on the \code{full} and \code{last} arguments.  
  Most calculations follow descriptions in Seber(1982), pages 414-418. If only two ages are present, a warning message 
  is generated and the catch curve method is not calculated. Plus groups are not allowed. Any NAs represent no estimates due to some issue with model fit
  like convergence. If age samples were collected via a survey using gears such as seines or trawl, or 
  were subsampled from catch, the least biased estimators are the "pois" and "crcb" methods (Nelson, 2019). 
}
\value{
  \item{results}{list element containing table of parameters and standard errors.}
  \item{data}{list element containing the age frequency data used in the analysis.}
}
\references{
Seber, G. A. F. 1982. The Estimation of Animal Abundance and Related Parameters, Second Edition. The Blackburn Press,
Caldwell, New Jersey. 654 pages. 

Maceina, M. J. and P. W. Bettoli. 1998. Variation in largemouth bass recruitment in four mainstream impoundments of the Tennessee River. 
N. Am. J. Fish. Manage. 18: 990-1003.

Millar, R. B. 2015. A better estimator of mortality rate from age-frequency data. Can. J. Fish. Aquat.
Sci. 72: 364-375.

Nelson, G. A. 2019. Bias in common catch-curve methods applied to age frequency data from fish surveys. ICES
J. Mar. Sci. doi:10.1093/icesjms/fsz085.

Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages.

Smith, M. W. and 5 others. 2012. Recommendations for catch-curve analysis. N. Am. J. Fish. Manage. 32: 956-967.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{data(rockbass)
agesurv(age=rockbass$age,full=6)}
\keyword{misc}
