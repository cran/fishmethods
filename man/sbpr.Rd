\name{sbpr}
\alias{sbpr}
\title{Spawning Stock Biomass-Per-Recruit Analysis}
\description{
  Spawning stock biomass-per-recruit(SBPR) analysis is conducted following Gabriel et al. (1989).  Reference points
of F and SBPR for a percentage of maximum spawning potential are calculated.}
\usage{
 sbpr(age = NULL, ssbwgt = NULL, partial = NULL, pmat = pmat,
 M = NULL, pF = NULL, pM = NULL, MSP = 40, plus = FALSE,
 oldest = NULL, maxF = 2, incrF = 1e-04, graph = TRUE)
}
\arguments{
  \item{age}{vector of cohort ages. If the last age is a plus group, do not add a "+" to the age.}
  \item{ssbwgt}{vector of spawning stock weights for each age. Length of vector must correspond to 
   the length of the age vector.}
  \item{partial}{partial recruitment vector applied to fishing mortality (F) to obtain partial F-at-age.  
    Length of this vector must match length of the age vector.}
  \item{pmat}{proportion of mature fish at each age. Length of this vector must match the length of the age vector.}
  \item{M}{vector containing a single natural mortality (M) rate if M is assumed constant over all ages, or a vector of
   Ms, one for each age. If the latter, the vector length match the length of the age vector.}
  \item{pF}{the proportion of fishing mortality that occurs before spawning.}
  \item{pM}{the proportion of natural mortality that occurs before spawning.}
  \item{MSP}{the percentage of maximum spawning potential (percent MSP reference point) for which F and SBPR should be 
    calculated.}
  \item{plus}{a logical indicating whether the last age is a plus-group. Default=FALSE.}
  \item{oldest}{if plus=TRUE, a numeric value indicating the oldest age in the plus group. }
  \item{maxF}{the maximum value of F range over which SBPR will be calculated. SBPR is calculated for F = 0 to maxF.}
  \item{incrF}{F increment for SBPR calculation.}
  \item{graph}{a logical indicating whether SPR and Percent Max SPR versus F should be plotted. Default=TRUE.}
}
\details{
  Spawning stock biomass-per-recruit analysis is conducted following Gabriel et al. (1989).  The F and SBPR for the 
  percentage maximum spawning potential reference point are calculated. If the last age is a plus-group, the cohort 
  is expanded to the \code{oldest} age and the \code{ssbwgt}, \code{partial}, \code{pmat}, and \code{M} values 
  for the plus age are applied to the expanded cohort ages.
}
\value{
  \item{Reference_Points}{F and SBPR values for the percentage MSP}
  \item{SBPR_vs_F }{Spawning stock biomass-per-recruit values for each F increment}
}
\references{Gabriel, W. L., M. P. Sissenwine, and W. J. Overholtz. 1989. Analysis of spawning stock biomass 
  per recruit: an example for Georges Bank haddock. North American Journal of Fisheries Management 9: 383-391.}

\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{ypr}}}
\examples{
data(haddock)
sbpr(age=haddock$age,ssbwgt=haddock$ssbwgt,partial=haddock$partial,
pmat=haddock$pmat,M=0.2,pF=0.2, pM=0.1667,MSP=30,plus=FALSE,maxF=2,
incrF=0.001)		
}
\keyword{misc}
