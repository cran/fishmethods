\name{ypr}
\alias{ypr}
\title{Yield-Per-Recruit Analysis}
\description{
 Yield-per-recruit (YPR) analysis is conducted following the modified Thompson-Bell algorithm. 
}
\usage{
ypr(age = NULL, wgt = NULL, partial = NULL, M = NULL,
 plus = FALSE, oldest = NULL, maxF = 2, incrF = 0.001, options = c(1), Fsol = NULL,
 graph = TRUE)
}
\arguments{
  \item{age}{the vector of cohort ages, e.g. c(1,2,3,4,5). If the last age is a plus group, 
do not add a "+" to the age.}
  \item{wgt}{the vector of catch weights for each age, e.g. c(0.2,0.4,0.7,1.0,1.2). Length of vector must correspond to the length of the age vector.}
  \item{partial}{the partial recruitment vector applied to fishing mortality (F) to obtain partial F-at-age.  Length of the partial recruitment vector
must correspond to the length of the age vector. }
  \item{M}{vector containing a single natural mortality (M) rate if M is assumed constant over all ages, or a vector of
 Ms, one for each age. If the latter, the vector length must correspond to the length of the age vector. }
  \item{plus}{a logical value indicating whether the last age is a plus-group. Default is FALSE.}
  \item{oldest}{if plus=TRUE, a numeric value indicating the oldest age in the plus group. }
  \item{maxF}{the maximum value of F range over which YPR will be calculated. YPR is calculated for F = 0 to maxF.}
  \item{incrF}{ F increment for YPR calculation.}
  \item{options}{1 = generate yield-per-recruit values for F ranging from 0 to maxF by incrF.
     2 =find a single YPR value for a given value of F (Fsol).Default = c(1).}
  \item{Fsol}{F for which to obtain a corresponding YPR. Default = NULL.}
  \item{graph}{logical indicating whether YPR versus F should be plotted. Default=TRUE.}
}
\details{
  Yield-per-recruit analysis is conducted following the modified Thompson-Bell algorithm.  Reference points
Fmax and F0.1 are calculated. If the last age is a plus-group, the cohort is expanded to the
\code{oldest} age and the \code{wgt}, \code{partial}, and \code{M} values for the plus age are applied to the expanded cohort ages.
}
\value{
  \item{Reference_Points}{If options = 1, F and yield-per-recruit values for Fmax and F0.1}
  \item{YPR_at_Fsol}{If options = 2, YPR at corresponding \code{Fsol}.}
  \item{F_vs_YPR}{For options = 1, yield-per-recruit values for each F increment.}

}
\references{Gabriel, W. L., M. P. Sissenwine, and W. J. Overholtz. 1989. Analysis of spawning stock biomass per recruit:
an example for Georges Bank haddock. North American Journal of Fisheries Management 9: 383-391.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{sbpr}}}
\examples{
data(haddock)
ypr(age=haddock$age,wgt=haddock$ssbwgt,partial=haddock$partial,M=0.4,
plus=TRUE,oldest=100,Fsol=0.2,maxF=2,incrF=0.01)
}
\keyword{misc}
