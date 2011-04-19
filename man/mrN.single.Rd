\name{mrN.single}
\alias{mrN.single}
\title{Estimate of Population Size from a Single Mark-Recapture Experiment}
\description{
 Estimates population sizes, standard errors, and confidence intervals for 
the bias-corrected Petersen and the Bailey binomial estimators.
}
\usage{
mrN.single(M = NULL, C = NULL, R = NULL, alpha = 0.05)
}
\arguments{
  \item{M}{ Number of marked animals released}
  \item{C}{Number of animals captured}
  \item{R}{ Number of animals recaptured}
  \item{alpha}{ alpha level for confidence intervals}
}
\details{
  The  bias-corrected Petersen estimator and its variance (Seber 2002: p.60),  and the Bailey binomial estimator and its variance 
(Seber 2002: p.61) are calculated.  The hypergeometric distribution is used to estimate confidence intervals for the Petersen model 
and the binomial distribution is used to estimate confidence intervals for the Bailey model.  
}
\value{Dataframe containing the population estimates (N), standard errors of N, the lower confidence limits (LCI),
and the upper confidence limits(UCI).
}
\references{ Seber, G. A. F. 2002. \emph{The Estimation of Animal Abundance and Related Parameters, Second Edition}. The Blackburn Press,
Caldwell, New Jersey. 654 p.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
mrN.single(M=948,C=421,R=167)
}
\keyword{misc}