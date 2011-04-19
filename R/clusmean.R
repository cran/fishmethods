################################################################################
#
#      Calculations of Simple Cluster Mean Length Using Package 'survey'
#          and Using Formulae in Pennington et. al (2002).
#
#        Example: All Fish in Tow Measured
################################################################################

#data<-read.csv("P:/Rwork/book/cod.csv",header=TRUE)

clusmean<-function(popchar=NULL, cluster=NULL, clustotal=NULL){
     means<-aggregate(popchar,list(cluster),mean)
        names(means)<-c("cluster","mean")
     sumM<-aggregate(popchar,list(cluster,clustotal),length)
        names(sumM)<-c("cluster","M","m")
     stats<-merge(means,sumM,by.x=c("cluster"),by.y=c("cluster"))
     n<-length(stats$cluster)
     R<-sum(stats$mean*stats$M)/sum(stats$M)
     varR<-sum(((stats$mean-R)^2*(stats$M/mean(stats$M))^2)/(n*(n-1)))
     rss<-aggregate((popchar-R)^2,list(cluster),sum)
        names(rss)<-c("cluster","rss")
     stats2<-merge(stats,rss,by.x="cluster",by.y="cluster")
     s2x<-sum((stats2$M/stats2$m)*stats2$rss)/(sum(stats2$M)-1)
     M<-sum(stats2$M);m<-sum(stats2$m);meff<-round(s2x/varR,0)
     ans<-matrix(NA,1L,7L)
     ans<-cbind(n,M,m,R,varR,s2x,meff)
     dimnames(ans)[[1]][1]<-list(" ")
     return(ans) 
}
  
#dd<-clusmean(popchar=data$len,cluster=data$station,clustotal=data$total)

