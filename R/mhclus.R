mhclus<-function(popchar=NULL, cluster=NULL, mhtype=c(1,2)){
     if(is.null(cluster)) stop("cluster vector not specified.")
     if(is.null(popchar)) stop("popchar vector not specified.")
     if(any(c(length(cluster),length(popchar)) %in% length(cluster)=="FALSE"))
         stop("vectors lengths are different.")
      if(!is.numeric(popchar)) stop("popchar must be numeric.")
       temp<-NULL
       temp<-as.data.frame(cbind(cluster,popchar))
       names(temp)<-c("cluster","popchar")
       temp<-temp[!is.na(temp$cluster) & !is.na(temp$popchar),]
       n<-length(mhtype)
       out<-matrix(NA,n,1L)
       dimnames(out)<-list(rep(NA,n),c("value"))
       cnt<-0
 if(any(mhtype==1)){
     # ICC
       cnt<-cnt+1
       dimnames(out)[[1]][cnt]<-list("Intra-cluster rho")
      ssto<-sum((temp$popchar-mean(temp$popchar))^2)
      myui<-aggregate(temp$popchar,list(temp$cluster),mean)
      names(myui)<-c("cluster","mean")
      newdats<-merge(temp,myui,by.x="cluster",by.y="cluster")
      ssw<-sum((newdats$popchar-newdats$mean)^2)
      catch<-aggregate(temp$popchar,list(temp$cluster),length)
      avgcatch<-mean(catch[,2])
      out[cnt,1]<-1-avgcatch/(avgcatch-1)*(ssw/ssto)
 }
if(any(mhtype==2)){
    # adjusted r2
     cnt<-cnt+1
     dimnames(out)[[1]][cnt]<-list("Adjusted r-square")
    mm<-aggregate(temp$popchar,list(temp$cluster),length)
    sm<-sum(mm[,2]-1)
    msw<-ssw/sm
    S2<-ssto/(sum(mm[,2])-1)
    out[cnt,1]<-1-msw/S2
 }
  return(out)
}
