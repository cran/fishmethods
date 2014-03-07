\name{clus.t.test}
\alias{clus.t.test}
\title{Correcting a Two-Sample Test for Clustering}
\description{
 Calculates Hedges (2007) t-statistic adjustment and degrees of freedom 
 for a t-test assuming unequal variances and clustered data with clusters 
 of unequal size.
}
\usage{clus.t.test(popchar = NULL, cluster = NULL, group = NULL,
      rho = NULL, alpha = 0.05, alternative = c("two.sided"))
}
\arguments{
  \item{popchar}{vector of population characteristic measurements (e.g., length, weight, etc.). 
    One row represents the measurement for an individual.}
  \item{cluster}{vector of numeric or character codes identifying individual clusters (or hauls).}
  \item{group}{vector of group membership identifiers.}
  \item{rho}{common intra-cluster correlation for groups.}
  \item{alpha}{alpha level used to calculate t critical value. Default=0.05}
  \item{alternative}{a character string specifying the alternative hypothesis, must be 
    one of "two.sided" (default), "greater" or "less".}
}
\details{
 A two-sample t-test with unequal variances (Sokal and Rohlf, 1995) is performed on clustered data. 
 The t-statistic and degrees of freedom are corrected for clustering according to 
 Hedges (2007). 
}
\value{List with null hypothesis of test and matrix table with mean of each group, rho, 
 ntilda (Equation 14 of Hedges 2007), nu (Equation 15), degrees of freedom (Equation 16), 
 uncorrected t-statistic, cu (Equation 18), the t-statistic adjusted for clustering, critical t value for 
 degrees of freedom and alpha, and probability of significance.
}
\references{Sokal,R.R.and F.J.Rohlf. 1995. Biometry, 3rd Edition. W.H. Freeman and Company, New York, 
            NY. 887 p.

  Hedges,L.V. 2007. Correcting a significance test for clustering. Journal of Educational and Behavioral 
  Statistics. 32: 151-179.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
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
 clus.t.test(popchar=newdata$length,cluster=newdata$station,
            group=newdata$region,rho=0.72,
            alpha=0.05,alternative="two.sided") 
}
\keyword{misc}
