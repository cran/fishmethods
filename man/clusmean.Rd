\name{clusmean}
\alias{clusmean}
\title{Estimation of Population Attributes and Effective Sample Size 
for Fishes Collected Via Cluster Sampling}
\description{
Calculates mean attribute, variance, and effective sample size for samples collected by simple random cluster sampling.
}
\usage{
clusmean(popchar = NULL, cluster = NULL, clustotal = NULL)
}
\arguments{
  \item{popchar}{vector of population characteristic measurements (e.g., length, weight, etc.). One row represents the measurement for an individual. }
  \item{cluster}{vector of numeric or character codes identifying individual clusters (or hauls).}
  \item{clustotal}{vector of total number of fish caught per cluster.}
}
\details{
In fisheries, gears (e.g., trawls, haul seines, gillnets, etc.) are used to collect fishes. Often, estimates of mean population attributes (e.g., mean length) are desired.
The samples of individual fish are not random samples, but cluster samples because the "haul" is the primary sampling unit.  Correct estimation of mean attributes requires the use of cluster sampling formulae.  
Estimation of the general mean attribute and variance follows Pennington et al. (2002). In addition, the effective sample size (the number of fish that would need to be sampled randomly to obtained the same precision 
as the mean estimate from cluster sampling) is also calculated.  The total number of fish caught in a cluster (\code{clustotal}) allows correct computation
for one- and two-stage sampling of individuals from each cluster (haul). 
}
\value{
 Matrix table of total number of clusters (n), total number of samples (M), total number of samples
measured (m), the mean attribute (R), variance of R (varR), variance of population attribute (s2x), and effective sample size (meff).
}
\references{ Pennington, M. , L. Burmeister, and V. Hjellvik. 2002. \emph{Assessing the precision of frequency distributions estimated from trawl-survey samples}
Fish. Bull. 100:74-80.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
data(codlengths)
clusmean(popchar=codlengths$tl,cluster=codlengths$station,
clustotal=codlengths$total)
}
\keyword{misc}
