\name{clus.rho}
\alias{clus.rho}
\title{Intracluster Correlation Coefficients for Clustered Data
}
\description{Calculates the intracluster correlation coefficients according to 
  Lohr (1999) and Donner (1986) for a single group
}
\usage{
 clus.rho(popchar=NULL, cluster = NULL, type = c(1,2,3), est = 0, nboot = 500)
}
\arguments{
  \item{popchar}{vector containing containing the population characteristic (e.g., length, weight, etc.). 
       One line per individual.}
  \item{cluster}{vector containing the variable used to identify the cluster. Identifier can be numeric or character.}
  \item{type}{method of intracluster correlation calculation. 1 = Equation 5.8 of Lohr (1999),
       2 = Equation 5.10 in Lohr (1999) and 3 = ANOVA.  Default = c(1,2,3).}
  \item{est}{estimate variance and percentiles of intracluster correlation coefficients via boostrapping.
       0 = No estimation (Default), 1 = Estimate.}
  \item{nboot}{number of boostrap replicates for estimation of variance. nboot = 500 (Default).}
}

\details{
     The intracluster correlation coefficient (rho) provides a measure of similarity within clusters.
     \emph{type} = 1  is defined to be the Pearson correlation coefficient for NM(M-1)pairs (yij,yik) for i 
     between 1 and N and j<>k (see Lohr (1999: p. 139). The average cluster size is used as the equal cluster 
     size quantity in Equation 5.8 of Lohr (1999). If the clusters are perfectly homogeneous (total variation is all 
     between-cluster variability), then ICC=1.

     \emph{type} = 2 is the adjusted r-square, an alternative quantity following Equation 5.10 in Lohr (1999). It is the
     relative amount of variability in the population explained by the cluster means, adjusted for the number 
     of degrees of freedom. If the clusters are homogeneous, then the cluster means are highly variable relative 
    to variation within clusters, and the r-square will be high.  
     
     \emph{type} = 3 is calculated using one-way random effects models (Donner, 1986).
     The formula is
         
     rho = (BMS-WMS)/(BMS+(m-1)*WMS)

     where BMS is the mean square between clusters, WMS is the mean square within clusters and m is the
     adjusted mean cluster size for clusters with unequal sample size. All clusters with zero elementary 
     units should be deleted before calculation. \emph{type} = 3 can be used with binary data 
    (Ridout et al. 1999)

     If \emph{est}=1, the boostrap mean (value), variance of the mean and 0.025 and 0.975 percentiles are outputted. 
}
\value{rho values, associated statistics, and bootstrap replicates}
\references{ Donner, A. 1986. A review of inference procedures for the intraclass correlation coefficient
 in the one-way random effects model. International Statistical Review. 54: 67-82.

 Lohr, S. L. Sampling: design and analysis. Duxbury Press,New York, NY. 494 p.

 Ridout, M. S., C. G. B. Demetrio, and D. Firth. 1999. Estimating intraclass correlation for
  binary data. Biometrics 55: 137-148.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{clus.lf} \link{clus.str.lf} \link{clus.mean}}}
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
 clus.rho(popchar=outs$len,cluster=outs$tow)}
\keyword{misc}
