\name{mhclus}
\alias{mhclus}
\title{Measures of Homogeneity Within Clusters
}
\description{Calculates measures of homogeneity within clusters according to Lohr (1999)
}

\usage{
 mhclus(popchar=NULL, cluster = NULL, mhtype=c(1,2))
}
\arguments{
  \item{popchar}{vector containing containing the population characteristic (e.g., length, weight, etc.). 
       One line per individual.}
  \item{cluster}{vector containing the variable used to identify the cluster. Identifier can be numeric or character.}
  \item{mhtype}{the type of homogeneity measure to calculate: 1 = intra-cluster correlation coefficient, 
    2 = adjusted r-square.Default = c(1,2).}
}

\details{
     The intracluster correlation coefficient (ICC) and adjusted r-square provides measures of homeogeneity within clusters.  
     ICC is defined to be the Pearson crrelation coefficient for NM(M-1)pairs (yij,yik) for i between 1 and N and j<>k (see 
     Lohr (1999: p. 139). The average cluster size is used as the equal cluster size quantity in Equation 5.8 of Lohr (1999).
     If the clusters are perfectly homogeneous (total variation is all between-cluster variability), then ICC=1. 

     Adjusted r-square is an alternative quantity following Equation 5.10 in Lohr (1999). It is the relative amount of variability
     in the population explained by the cluster means, adjusted for the number of degrees of freedom.  If the clusters are homogeneous,
     then the cluster means are highly variable relative to variation within clusters, and the r-square will be high.
 
     All clusters with zero elementary units should be deleted before calculation.
  
}
\value{A matrix of values by mhtype.
}
\references{ Lohr, S. L. Sampling: design and analysis. Duxbury Press,New York, NY. 494 p.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{lfstrclus} \link{lfclus}}}
\examples{
  data(codcluslen)
  tem<-codcluslen[codcluslen[,1]=="NorthCape" & codcluslen[,3]>0,]
  outs<-data.frame(tow=NA,len=NA)
  cnt<-0
  for(i in 1:as.numeric(length(tem$number))){
    for(j in 1:tem$number[i]){
     cnt<-cnt+1
     outs[cnt,1]<-tem$tow[i]
     outs[cnt,2]<-tem$length[i]
   }
 }
 mhclus(popchar=outs$len,cluster=outs$tow)}
\keyword{misc}
