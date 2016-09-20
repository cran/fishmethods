\name{catch.select}
\alias{catch.select}
\title{Selectivity Ogive from a Catch Curve}
\description{
  Estimates selectivity-at-length from catch lengths and von Bertalanffy growth parameters.
}
\usage{
catch.select(len = NULL, lenmin = NULL, binsize = NULL, 
peakplus = 1, Linf = NULL, K = NULL, t0 = NULL, subobs = FALSE)
}
\arguments{
  \item{len}{vector of lengths. One row per individual.}
  \item{lenmin}{the starting length from which to construct length intervals.}
  \item{binsize}{the length interval width.  Must be >0. This is used to create the lower length of intervals starting from \code{lenmin} to the maximum observed in \code{len}. 
}
  \item{peakplus}{numeric. Allows user to specify the number of length intervals following the length interval at the peak log(catch/deltat) to use as the start length interval in the catch curve analysis. Default = 1 based on standard catch curve analysis recommendations of Smith et al. (2012).
}
 \item{Linf}{numeric. The L-infinity value from a von Bertalanffy growth equation. This is a required value.}
  \item{K}{numeric. The growth coefficient from a von Bertalanffy growth equation. This is a required value.}
  \item{t0}{numeric. The t-subzero value from a von Bertalanffy growth equation. This is a required value.}
  \item{subobs}{logical. If the "observed" selectivity for those under-represented length intervals not used in the catch curve analysis is equal to 1, the inverse logit (used in fit of selectivity ogive) can not be calculated. If \code{subobs} is set to TRUE, 1 will be substituted with 0.9999}
}
\details{
 This function applies the method of Pauly (1984) for calculating the selectivity-at-length from catch lengths and parameters from a von Bertalanffy growth curve. See Sparre and Venema(1998) for a detailed example of the application.

Length intervals are constructed based on the \code{lenmin} and \code{binsize} specified, and the maximum length observed in the data vector. Catch-at-length is tabularized using the lower and upper intervals and the data vector of lengths. The inclusion of a length in an interval is determined by lower interval>=length<upper interval. The age corresponding to the interval midpoint (\code{t}) is determined using the von Bertalanffy equation applied to the lower and upper bounds of each interval, summing the ages and dividing by 2. \code{deltat} is calculated for each interval using the equation: (1/k)*log((Linf-L1)/(Linf-L2)) where L1 and L2 are the lower and upper bounds of the length interval. \code{log(catch/deltat)} is the dependent variable and \code{t} is the predictor used in linear regression to estimate Z. Using the parameters from the catch curve analysis,  "observed" selectivities (\code{stobs}) for the length intervals not included in the catch curve analysis are calculated using the equation: stobs=catch/(deltat*exp(a-Z*t)) where a and Z are the intercept and slope from the linear regression. The \code{stobs} values are transformed using an inverse logit (log(1/stobs-1)) and are regressed against \code{t} to obtain parameter estimates for the selectivity ogive.  The estimated selectivity ogive (\code{stest}) is then calculated as
stest=1/(1+exp(T1-T2*t)) where T1 and T2 are the intercept and slope from the log(1/stobs-1) versus \code{t} regression.

}
\value{list containing a dataframe with the lower and upper length intervals, the mid-point length interval, age corresponding to the interval mid-point, catch of the length interval, log(catch/deltat), the predicted log(catch/deltat) from the catch curve model fit (only for the peakplus interval and greater), the observed selectivities and the estimated selectivity, and two dataframes containing the parameters and their standard errors from the linear regressions for catch curve analysis and the selectivity ogive. 
}
\references{

Pauly, D. 1984.  Length-converted catch curves. A powerful tool for fisheries research in the tropics (Part III). ICLARM Fishbyte 2(1): 17-19.

Smith, M. W. and 5 others. 2012. Recommendations for catch-curve analysis. N. Am. J. Fish. Manage. 32: 956-967.

Sparre, P. and S. C. Venema. 1998. Introduction to tropical fish stock assessment. Part 1. Manual. FAO Fisheries Technical Paper, No. 206.1, Rev. 2. Rome. 407 p. Available on the world-wide web.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
data(sblen)
catch.select(len=sblen$len_inches,binsize=2,lenmin=10,peakplus=1,Linf=47.5,K=0.15, 
t0=-0.3)
}
\keyword{misc}

