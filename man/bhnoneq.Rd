\name{bhnoneq}
\alias{bhnoneq}
\title{Length-based Beverton-Holt Nonequilibrium Z Estimator}
\description{
  A nonequilibrium Beverton-Holt estimator of instantaneous total mortality (Z) from 
  length data.
}
\usage{
bhnoneq(year=NULL,mlen=NULL, ss=NULL, K = NULL, Linf = NULL, 
Lc = NULL, nbreaks = NULL, styrs = NULL, stZ = NULL, 
stsigma = NULL)
}
\arguments{
 \item{year}{the vector of year values associated with mean length data. The number 
             of year values must correspond to the number of length records. Include year value even if mean length and numbers (see below)
             are missing.}
  \item{mlen}{the vector of mean lengths for lengths >=Lc. One record for each year.}
 \item{ss}{the vector of numbers of observations associated with the mean length.}  
  \item{K}{the growth coefficient from a von Bertalanffy growth model.}
  \item{Linf}{the L-infinity coefficient from a von Bertalanffy growth model.}
  \item{Lc}{the length at first capture.}
  \item{nbreaks}{the number of times (breaks) mortality is thought to change over the time series. Must be 1 or greater}
  \item{styrs}{the starting guess(es) of the year(s) during which mortality is thought to change. 
       The number of starting guesses must match the number of mortality breaks, should be separated by
       commas within the concatentation function and should be within the range of years present 
       in the data. 
        }
  \item{stZ}{the starting guesses of Z values enclosed within the concatentation function. There should be \emph{nbreaks+1} values provided. 
        }
\item{stsigma}{the starting guess of sigma.} 
}
\details{
  The mean lengths for each year for lengths>=Lc. Following Gedamke and Hoening(2006), the model estimates \code{nbreaks+1} Z values, the year(s) in which the
  changes in mortality began, the standard deviation of lengths>=Lc, and standard errors of all parameters. An AIC value is produced for model comparison.
  The estimated parameters for the number of \code{nbreaks} is equal to \code{2*nbreaks+2}.  Problematic parameter estimates may have extremely large t-values or
  extremely small standard error.  Try different starting values to ensure consistent parameter estimates.
}
\value{
  \item{results}{list element containing table of parameters with estimates, 
   standard errors, and t-values.}
  \item{obs}{list element containing year, mean length of lengths >=Lc, 
   and sample size.}
  \item{pred}{list element containing year and predicted mean lengths.}
 }
\references{Gedamke, T. and J. M. Hoenig. 2006. Estimating mortality from mean length
  data in nonequilibrium situations, with application to the assessment of goosefish. Trans. Am. Fish. Soc. 135:476-487}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\note{ Todd Gedamke provided the predicted mean length code in C++.
}
\seealso{\code{\link{bheq2}} \code{\link{bheq2}}}
\examples{data(goosefish)
bhnoneq(year=goosefish$year,mlen=goosefish$mlen, ss=goosefish$ss,
K=0.108,Linf=126,Lc=30,nbreaks=1,styrs=c(1982),stZ=c(0.1,0.3),
stsigma=20)
}
\keyword{misc}
