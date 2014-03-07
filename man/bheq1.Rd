\name{bheq1}
\alias{bheq1}
\title{Length-based Beverton-Holt Equilibrium Total Instantaneous Mortality Estimator}
\description{
  The equilibrium Beverton-Holt estimator of instantaneous total mortality (Z) from 
  length data with bootstrapped standard errors
}
\usage{
bheq1(len, K = NULL, Linf = NULL, Lc = NULL, nboot = 100)
}
\arguments{
  \item{len}{the vector of length data. Each row represents one record per individual fish.}
  \item{K}{the growth coefficient from a von Bertalanffy growth model.}
  \item{Linf}{the L-infinity coefficient from a von Bertalanffy growth model.}
  \item{Lc}{the length at first capture.}
  \item{nboot}{the number of bootstrap runs. Default=100.}
}
\details{
  The standard Beverton-Holt equilibrium estimator of instantaneous total mortality (Z)
  from length data (page 365 in Quinn and Deriso (1999)) is calculated. The mean length 
  for lengths >=Lc is calculated automatically. Missing data are removed prior to calculation.
 Estimates of standard error are made by bootstrapping length data >=Lc using package \code{boot}. 
}
\value{Dataframe of length 1 containing mean length>=Lc, sample size>=Lc, Z estimate and standard error.
}
\references{Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{bheq2}} \code{\link{bhnoneq}}}
\examples{data(pinfish)
bheq1(pinfish$sl,K=0.33,Linf=219.9,Lc=120,nboot=200)
}
\keyword{misc}
