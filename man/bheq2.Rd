\name{bheq2}
\alias{bheq2}
\title{Length-based Beverton-Holt Equilibrium Total Instantaneous Mortality Estimator with Bias-Correction}
\description{
  The equilibrium Beverton-Holt estimator of instantaneous total mortality (Z) from length data using Ehrhardt and Ault (1992) bias-correction
}
\usage{
bheq2(len = NULL, Linf = NULL, K = NULL, Lc = NULL, La = NULL,
 nboot = 100)
}
\arguments{
  \item{len}{the vector of length data. Each row represents one record per individual fish.}
  \item{K}{the growth coefficient from a von Bertalanffy growth model.}
  \item{Linf}{the L-infinity coefficient from a von Bertalanffy growth model.}
  \item{Lc}{the length at first capture.}
  \item{La}{the largest length of the largest size class.}
  \item{nboot}{the number of bootstrap runs. Default=100.}
}
\details{
  The  Beverton-Holt equilibrium estimator of instantaneous total mortality (Z)
  from length data (page 365 in Quinn and Deriso (1999)) using Ehrhardt and Ault (1992) bias-correction is calculated. The mean length 
  for lengths >=Lc is calculated automatically. Missing data are removed prior to calculation.
 Estimates of standard error are made by bootstrapping length data >=Lc using package \code{boot}. 
}
\value{Dataframe of length 1 containing mean length>=Lc, sample size>=Lc, Z estimate and standard error.
}
\references{Quinn, T. J. and R. B. Deriso. 1999. Quantitative Fish Dynamics. Oxford University Press, New York, New York. 542 pages

Ehrhardt, N. M. and J. S. Ault. 1992. Analysis of two length-based mortality models applied to bounded
catch length frequencies. Trans. Am. Fish. Soc. 121:115-122.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{bheq1}} \code{\link{bhnoneq}}}
\examples{data(pinfish)
bheq2(pinfish$sl,K=0.33,Linf=219.9,Lc=120,La=255,nboot=200)
}
\keyword{misc}

