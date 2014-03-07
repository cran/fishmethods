\name{clus.mean}
\alias{clus.mean}
\title{Estimation of Population Attributes and Effective Sample Size 
for Fishes Collected Via Cluster Sampling}
\description{
Calculates mean attribute, variance, effective sample size, and degrees of freedom for samples collected 
 by simple random cluster sampling.
}
\usage{
clus.mean(popchar = NULL, cluster = NULL, clustotal = NULL, rho = NULL, nboot = 1000)
}
\arguments{
  \item{popchar}{vector of population characteristic measurements (e.g., length, weight, etc.). 
    One row represents the measurement for an individual. }
  \item{cluster}{vector of numeric or character codes identifying individual clusters (or hauls).}
  \item{clustotal}{vector of total number of fish caught per cluster.}
  \item{rho}{intracluster correlation coefficient for data. If NULL, degrees of freedom are not calculated.}
  \item{nboot}{number of bootstrap samples for calculation of bootstrap variance. Default = 1000}
}
\details{
In fisheries, gears (e.g., trawls, haul seines, gillnets, etc.) are used to collect fishes. Often, estimates of mean population attributes (e.g., mean length) are desired.
The samples of individual fish are not random samples, but cluster samples because the "haul" is the primary sampling unit.  Correct estimation of mean attributes requires 
the use of cluster sampling formulae.  Estimation of the general mean attribute and usual variance approximation follows Pennington et al. (2002). 
Variance of the mean is also estimated using the jackknife and bootstrap methods (Pennington and Volstad, 1994; Pennington et al., 2002).
In addition, the effective sample size (the number of fish that would need to be sampled randomly to obtained the same precision 
as the mean estimate from cluster sampling) is also calculated for the three variance estimates.  The total number of fish caught in a cluster 
(\code{clustotal}) allows correct computation for one- and two-stage sampling of individuals from each cluster (haul).
In addition, if rho is specified, degrees of freedom are calculated by using Hedges (2007) for unequal cluster sizes (p. 166-167). 
}
\value{
 Matrix table of total number of clusters (n), total number of samples (M), total number of samples
measured (m), the mean attribute (R), usual variance approximation (varU), jackknife variance (varJ), bootstrap variance (varB),
variance of population attribute (s2x), usual variance effective sample size (meffU), jackknife variance effective sample size,
(meffJ), bootstrap variance effective sample size (meffB) and degrees of freedom (df) if applicable. 
}
\references{Hedges,L.V. 2007. Correcting a significance test for clustering. Journal of Educational and Behavioral 
  Statistics. 32: 151-179.

Pennington, M. and J. J. Volstad. 1994. Assessing the effect of intra-haul correlation and variable density 
on estimates of population characteristics from marine surveys.  Biometrics 50: 725-732.

Pennington, M. , L. Burmeister, and V. Hjellvik. 2002. Assessing the precision of frequency distributions estimated 
from trawl-survey samples. Fish. Bull. 100:74-80.

}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
data(codcluslen)
temp<-codcluslen[codcluslen$region=="NorthCape" & codcluslen$number>0,]
temp$station<-c(paste(temp$region,"-",temp$tow,sep=""))
total<-aggregate(temp$number,list(temp$station),sum)
names(total)<-c("station","total")
temp<-merge(temp,total,by.x="station",by.y="station")
newdata<-data.frame(NULL)
for(i in 1:as.numeric(length(temp[,1]))){
  for(j in 1:temp$number[i]){
    newdata<-rbind(newdata,temp[i,])
  }
}
clus.mean(popchar=newdata$length,cluster=newdata$station,
         clustotal=newdata$total)
}
\keyword{misc}
