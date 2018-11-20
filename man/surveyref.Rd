\name{surveyref}
\alias{surveyref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantitative reference points from stock abundance indices based on research surveys}
\description{
 This function implements the methodology of Helser and Hayes (1995) for generating quantitative reference points from relative 
abundance indices based on research surveys
}
\usage{
surveyref(x = NULL, refpt = 25, compyear = NULL, reffix = FALSE,
 refrange = NULL, nboot = 500, allboots = FALSE, nreps = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{output object from function \code{surveyfit}.}
  \item{refpt}{the lower quantile (percentile) of the fitted time series used as the reference point.}
  \item{compyear}{ the index year to compare to the reference point.  Multiple years can be included in the comparison using the \code{c()} function.}
  \item{reffix}{a logical value specifying whether the lower quantile should be determined from a fixed set of years.
 Default = FALSE.}
  \item{refrange}{If \code{reffix} = TRUE, the beginning and ending year of the time series to include in 
determination of the lower quantile.  The values should be enclosed within \code{c()} (e.g., c(1963,1983)).}
  \item{nboot}{the number of bootstrap replicates.}
  \item{allboots}{a logical value specifying whether the fitted values for the bootstrap replicates 
should be included in the output. Default = FALSE.}
  \item{nreps}{the number of samples to draw  in function \code{pgen}. Default = 10000.}
}
\details{
  Using the output object from function \code{surveyfit}, the methodology of Helser and Hayes (1995) is applied to 
generate the probability distribution that the abundance index value for a given year lies below the value of a 
lower quantile (reference point).  The procedure is : 1) add to the original fitted time series residuals randomly selected
 with replacement from the Pennington model fit, 2) repeat this \code{nboot} times to create new time series, 
3) fit the Pennington model to each new time series using the original theta estimate to get \code{nboot} replicates
 of new fitted time series, and 4) determine the lower quantile for each new fitted time 
series.  The probability of the abundance index being less than the quartile reference point is calculated using 
function \code{pgen} with \code{comp}=1. 

If comparisons between the current year's index and the reference point will be made year-after-year, Helser and Hayes
 (1995) recommend using a fixed set of years to select the lower quantile.  This procedure will avoid a change in 
reference point over time as a survey time series is updated. Use arguments \code{reffix} and \code{refrange} to 
accomplish this.
}
\value{list containing the lower quantile of the original fitted time series and the mean quantile of the
fitted bootstrap replicates (\code{comp_refpt}),  the original fitted time series values versus the mean of the fitted
 bootstrap time series values(\code{comp_fitted}), the empirical distribution of the selected index (\code{emp_dist_index}),
 the empirical distribution of the lower quantile (\code{emp_dist_refpt}), the probability that the index
 value lies below the reference point for a given decision confidence level (\code{prob_index}), and, if argument \code{allboots} is TRUE, the fitted values
 of the bootstrap replicates (\code{boot_runs}).
}
\references{Helser, T. E. and D. B. Hayes. 1995.  Providing quantitative management advice from stock abundance 
indices based on research surveys. Fishery Bulletin 93: 290-298.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{surveyfit}}} 
\examples{
data(wolffish)
out<-surveyfit(year=wolffish$year,index=wolffish$index,logtrans=TRUE)
surveyref(out,refpt=25,compyear=c(1990))
}
\keyword{misc}

