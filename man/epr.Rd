\name{epr}
\alias{epr}
\title{Eggs-Per-Recruit Analysis}
\description{
  Eggs-per-recruit(EPR) analysis is conducted following Gabriel et al. (1989) except fecundity-at-age is substituted for weight-at-age.  Reference points
of F and EPR for percentage of maximum spawning potential are calculated.}
\usage{
epr(age = NULL, fecund = NULL, partial = NULL, pmat = pmat,
 M = NULL, pF = NULL, pM = NULL, MSP = 40, plus = FALSE,
 oldest = NULL, maxF = 2, incrF = 1e-04)
}
\arguments{
  \item{age}{vector of cohort ages. If the last age is a plus group, do not add a "+" to the age.}
  \item{fecund}{vector of fecundity (number of eggs per individual) for each age. Length of vector must correspond to the length of the age vector.}
  \item{partial}{partial recruitment vector applied to fishing mortality (F) to obtain partial F-at-age.  Length of this vector must match length of the age vector. }
  \item{pmat}{proportion of mature fish at each age. Length of this vector must match the length of the age vector.  }
  \item{M}{vector containing a single natural mortality (M) rate if M is assumed constant over all ages, or a vector of
  Ms, one for each age. If the latter, the vector length match the length of the age vector. }
  \item{pF}{the proportion of fishing mortality that occurs before spawning.}
  \item{pM}{the proportion of natural mortality that occurs before spawning.}
  \item{MSP}{the percentage of maximum spawning potential (percent MSP reference point) for which F and SBPR should be calculated.  }
  \item{plus}{a logical value indicating whether the last age is a plus-group. Default is FALSE.}
  \item{oldest}{if plus=TRUE, a numeric value indicating the oldest age in the plus group. }
  \item{maxF}{the maximum value of F range over which SBPR will be calculated. SBPR is calculated for F = 0 to maxF.}
  \item{incrF}{F increment for SBPR calculation.}
}
\details{
  Eggs-per-recruit analysis is conducted following Gabriel et al. (1989).  The F and EPR for the percentage maximum spawning potential reference point
  are calculated. If the last age is a plus-group, the cohort is expanded to the
  \code{oldest} age and the \code{fecund}, \code{partial}, \code{pmat}, and \code{M} values for the plus age are applied to the expanded cohort ages.
}
\value{
  \item{Reference_Points}{F and EPR values for the percentage MSP}
  \item{EPR_vs_F }{Eggs-per-recruit values for each F increment}
}
\references{Gabriel, W. L., M. P. Sissenwine, and W. J. Overholtz. 1989. Analysis of spawning stock biomass per recruit:
an example for Georges Bank haddock. North American Journal of Fisheries Management 9: 383-391.}

\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{ypr}} \code{\link{sbpr}}}
\examples{
    data(menhaden)
	epr(age=menhaden$age,fecund=menhaden$fecundity,partial=menhaden$partial,
	pmat=menhaden$pmat,M=menhaden$M,pF=0,pM=0,MSP=40,plus=TRUE,maxF=4,incrF=0.01,oldest=10)
}
\keyword{misc}
