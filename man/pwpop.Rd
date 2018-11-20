\name{pwpop}
\alias{pwpop}
\title{ Estimate Net Reproductive Rates Over Multiple Periods Of An Abundance Time Series Using Piecewise Regression}
\description{
  Function estimates net reproductive rates for periods of change over a time series of abundance data.
}
\usage{
pwpop(abund = NULL, year = NULL, periods = NULL, Cs = NULL,
 startR = NULL, upperR = NULL, lowerR = NULL, graph = TRUE)
}
\arguments{
  \item{abund}{the vector of time series of abundance data (e.g. run counts, indices of relative abundance, etc.).}
  \item{year}{the vector of years associated with abundance data.}
  \item{periods}{the number of periods over which to fit the population model.}
  \item{Cs}{the vector of user-specified initial starting value for year(s) of change - number of values equals \emph{periods} - 1 (enclose within c()).}
  \item{startR}{the vector of user-specified initial starting values for R - one value for each period (enclose within c()).}
  \item{upperR}{the vector of user-specified upper limits for R (one for each period) used in optimization (enclose within c()).}
  \item{lowerR}{the vector of user-specified lower limits for R (one for each period) used in optimization (enclose within c()).}
  \item{graph}{Logical specifying whether a graph of observed versus predicted values is plotted. Default=TRUE.}
}
\details{
  A simple population model is fitted to abundance data to estimate the net reproductive rate for specified periods of time.
  The model is Nt=N0*R^t where Nt is the abundance at time t, N0 is the estimated initial population size and R is the net
  reproductive rate.  R can be used as an indication that the population is stable (R=1), is increasing (R>1) or is declining
  (R<1) over a specified time period.  The fitted equation is the linearized form: \emph{log(Nt)=log(N0)+log(R)*t}, where log is
   the natural-log; therefore, zeros are not allowed.  
  
   To simultaneously estimate the parameters for periods of trends in the abundance data, a piecewise regression approach is
   used.  The linearized model is fitted separately to data for each period but models are linked so that the ending
   year for the preceding period is also the intercept for the current period. As an example, the models for three 
   periods are
    
    log(N1,t)=log(N1,0)+log(R1)*t for t<C1

    log(N2,t)=log(N1,0)+C1*(log(R1)-log(R2))+log(R2)*t for t>=C1 and t<C2

    log(N3,t)=log(N1,0)+C1*(log(R1)-log(R2))+C2*(log(R2)-log(R3))+log(R3)*t for t>=C2

   The parameters estimated for these models are log(N1,0), log(R1), C1, log(R2), C2, and log(R3). \emph{t} is time starting at 1 
   for the first year of abundance and ending at x for the last year of abundance(year information is still needed for 
   plotting). Entered Cs value are converted to the same scale as t. Back-transform the log(R) values using \emph{exp} 
   to obtain the R values for each period. The function \code{optim} is used to obtain parameter estimates and associated 
   standard errors by minimizing the sum of squares (log(N)-log(pred))^2. Add first year-1 to each C to put estimates on year scale.
}
\value{
  \item{Estimates}{list element with the parameter estimates and associated standard errors, residual sum of squares, Akaike's
   Information Criterion for least squares (AIC), and coefficient of determination (r2).}
  \item{Data}{list element with the abundance data, years, t, log predicted values, and back-transformation predicted values.}
}
\references{Neter, J. , M. H. Kutner, C. J. Nachtsheim, and W. Wasserman.  1996.  Applied Linear Statistical Models.
  The Magraw-Hill Companies. 1408 p.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
data(counts)
pwpop(abund = counts$number, year = counts$year,periods = 3, Cs = c(2000,2005), 
startR = c(0.5,0.5,0.5), 
upperR = c(10,10,10), 
lowerR = c(-10,-10,-10))
}
\keyword{misc}