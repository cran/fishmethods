\name{clus.rho.g}
\alias{clus.rho.g}
\title{Calculate A Common Intracluster Correlation Coefficient Among Groups
}
\description{Calculates a common intracluster correlation coefficients according to 
  Donner (1986: 77-79) for two or more groups with unequal cluster sizes, and tests for 
  homogeneity of residual error among groups and a common coefficient among groups.
}
\usage{
 clus.rho.g(popchar=NULL, cluster = NULL, group = NULL)
}
\arguments{
  \item{popchar}{vector containing containing the population characteristic (e.g., length, weight, etc.). 
       One line per individual.}
  \item{cluster}{vector containing the variable used to identify the cluster. Identifier can be numeric or character.}
  \item{group}{vector containing the identifier used for group membership of length data.  This variable is
  used to determine the number of groups.  Identifier can be numeric or character.}
}
\details{
     The intracluster correlation coefficient (rho) provides a measure of similarity within clusters.
     rho is calculated using a one-way nested random effects model (Donner, 1986: 77-79).
     The formula is
         
     rho = (BMS-WMS)/(BMS+(m-1)*WMS)

     where BMS is the mean square among clusters within groups, WMS is the mean square within clusters and m 
     is the adjusted mean cluster size for clusters with unequal sample sizes. All clusters with zero 
     elementary units should be deleted before calculation. In addition, approximate 95% confidence intervals 
     are generated and a significance test is performed.

    Bartlett's test is used to determine if mean square errors are constant among groups. If Bartlett's test is
    not significant, the test for a common correlation coefficient among groups is valid. 
    
}
\value{rho value and associate statistics}
\references{Donner, A. 1986. A review of inference procedures for the intraclass correlation coefficient
 in the one-way random effects model. International Statistical Review. 54: 67-82.
 }
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{clus.str.lf} \link{clus.lf} \link{clus.mean}}}
\examples{
  data(codcluslen)
   temp<-codcluslen[codcluslen$number>0,]
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
  newdata<-newdata[,-c(5)]
 clus.rho.g(popchar=newdata$length,cluster=newdata$station,group=newdata$region)
}
\keyword{misc}
