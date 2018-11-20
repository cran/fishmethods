\name{surveyfit}
\alias{surveyfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimating the Relative Abundance of Fish From a Trawl Survey}
\description{
 This function applies the time series method of Pennington (1986) for estimating relative abundance
 to a survey series of catch per tow data 
}
\usage{
surveyfit(year = NULL, index = NULL, logtrans = TRUE, graph = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{ vector containing the time series of numeric year labels.}
  \item{index}{vector containing the time series of mean catch per tow data.}
  \item{logtrans}{a logical value indicating whether the natural log-transform should be applied to the mean catch per tow values.
       Default is TRUE.}
 \item{graph}{a logical value indicating whether a graph of the observed and model fit
             should be drawn. Default is TRUE.}

}
\details{
  Parameters for a first difference, moving average model of order 1 are estimated from the trawl time series using function \code{arima}.
Following Equation 4 in Pennington (1986), fitted values are calculated from the model residuals and the estimate of theta.   
}
\value{List containing summary statistics (sample size (n), the first three sample autocorrelations (r1-r3) for the first differenced logged series)
 and parameter estimates (theta, theta standard error, and sigma2), the observed log-transformed index and fitted values, and the ARIMA function output.}
\references{Pennington, M. P. 1986. Some statistical techniques for estimating abundance indices from trawl surveys. Fishery Bulletin 84(3): 519-525.}
\author{ Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{surveyref}}} 
\examples{
data(yellowtail)
surveyfit(year=yellowtail$year,index=yellowtail$index)
}
\keyword{misc}

