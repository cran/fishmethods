\name{sbpr}
\alias{sbpr}
\title{Spawning Stock Biomass-Per-Recruit Analysis}
\description{
  Spawning stock biomass-per-recruit analysis (SSBPR) is conducted following Gabriel et al. (1989).}
\usage{
 sbpr(age=NULL,ssbwgt=NULL,partial=NULL,pmat=NULL,M=NULL,pF=NULL,pM=NULL,
               plus=FALSE,oldest=NULL,maxF=2,incrF=0.0001,options=c(1,2,3,4),MSP=NULL,
               SSBPR=NULL,Fsol=NULL,graph=TRUE)
}
\arguments{
  \item{age}{a numeric vector of cohort ages. If the last age is a plus group, do not add a "+" to the age.}
  \item{ssbwgt}{vector of spawning stock weights for each age. Length of vector must correspond to 
   the length of the age vector.}
  \item{partial}{partial recruitment vector applied to fishing mortality (F) to obtain partial F-at-age.  
    Length of this vector must match length of the age vector.}
  \item{pmat}{proportion of mature fish at each age. Length of this vector must match the length of the age vector.}
  \item{M}{vector containing a single natural mortality (M) rate if M is assumed constant over all ages, or a vector of
   Ms, one for each age. If the latter, the vector length match the length of the age vector.}
  \item{pF}{the proportion of fishing mortality that occurs before spawning.}
  \item{pM}{the proportion of natural mortality that occurs before spawning.}
  \item{plus}{a logical indicating whether the last age is a plus-group. Default=FALSE.}
  \item{oldest}{if plus=TRUE, a numeric value indicating the oldest age in the plus group. }
  \item{maxF}{the maximum value of F range over which SSBPR will be calculated. SSBPR is calculated for F = 0 to maxF.}
  \item{incrF}{F increment for SSBPR calculation.}
 \item{options}{1 = generate spawning stock biomass-per-recruit values for F ranging from 0 to maxF by incrF. 
        2 = find a single SSBPR value for a given value of F (Fsol), 3 = find F at a specified percent maximum SSBPR (MSP), 
        4 = find F for a given value of SSBPR. Default = c(1,2,3,4).}
  \item{Fsol}{F for which to obtain a corresponding spawning biomass-per-recruits value. Default = NULL.}
  \item{MSP}{the percentage of maximum spawning potential (percent MSP reference point) for which F and SBPR should be 
    determined. Default = NULL.}
  \item{SSBPR}{A spawning biomass-per-recruit value for which to obtain a corresponding F. Default = NULL.}
   \item{graph}{a logical indicating whether SSBPR versus F should be plotted for options = 1. Default=TRUE.}
}
\details{ 
 Spawning stock biomass-per-recruit analysis is conducted following Gabriel et al. (1989). If the last age is a 
  plus-group, the cohort is expanded to the \code{oldest} age and the \code{ssbwgt}, \code{partial}, \code{pmat},
  and \code{M} values for the plus age are applied to the expanded cohort ages. Multiple options are available
  to abbreviate calculations.
}
\value{
 \item{SSBPR_vs_F}{For option = 1, spawning stock biomass-per-recruit values for each F increment.}
  \item{SSBPR_at_Fsol}{If option = 2, the SSBPR value corresponding to Fsol.}
  \item{F_at_MSP}{If option = 3, the F reference point corresponding to MSP.}
  \item{F_at_SSBPR}{If option = 4, the F corresponding to a given SSBPR.}
}
\references{Gabriel, W. L., M. P. Sissenwine, and W. J. Overholtz. 1989. Analysis of spawning stock biomass 
  per recruit: an example for Georges Bank haddock. North American Journal of Fisheries Management 9: 383-391.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{ypr}}}
\examples{
data(haddock)
#Generate SSBPR versus F, and F for MSP reference point
sbpr(age=haddock$age,ssbwgt=haddock$ssbwgt,partial=haddock$partial,
pmat=haddock$pmat,M=0.2,pF=0.2, pM=0.1667,plus=FALSE,maxF=2,
incrF=0.001,MSP=30,options = c(1,3))		
}
\keyword{misc}
