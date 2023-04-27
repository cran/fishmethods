compare.lrt.plus<-function(...){
  foo.call <- as.character(match.call())[-1]
  all.x <- list(...)
  n <- length(all.x)
  if(n<=1) stop("Need at least two models to compare")
  names(all.x) <- unique(foo.call)
  check <- NULL
  for (i in 1:n){
    if(length(all.x[[i]])<6) stop(paste("Model ",i," fit failed"))
    if(!any(names(all.x[[i]])=="type")) stop("Only models from growthlrt.plus can be used in this function")
    if(all.x[[i]]$type!="growthlrt_plus") stop("Only models from growthlrt.plus can be used in this function")
  }  
  outpt<-data.frame(Model=1:n, df=rep(NA,n),AIC=rep(NA,n),BIC=rep(NA,n),logLik=rep(NA,n),Test=rep(NA,n),L.Ratio=rep(NA,n),p_value_LR=rep(NA,n),F.Stat=rep(NA,n),p_value_F=rep(NA,n))
 
  row.names(outpt) <- c(unique(foo.call))
  cnt <- 0
  for (x in list(...)) {
    cnt <- cnt + 1
    if(cnt==1) nobs<-nrow(x$residuals)
    outpt[cnt, 1] <- cnt
    outpt[cnt, 2] <- x$model_comp_df
    outpt[cnt, 3] <- x$results[row.names(x$results)=="AIC",1]
    outpt[cnt, 4] <- x$results[row.names(x$results)=="BIC",1]
    outpt[cnt, 5] <- x$results[row.names(x$results)=="LogLik",1]
  }
  
  for(i in 1:n) {
    if(i==1) outpt[i,6:10]<-""
    if(i>1){
      df<-abs(as.numeric(outpt[i-1,2])-as.numeric(outpt[i,2]))
      if(df<1) outpt[i,6:10]<-""
      if(df>=1){
        outpt[i,6]<-paste(i-1," vs ",i,sep="")
        outpt[i,7]<-round(2*abs((as.numeric(outpt[i-1,5])-as.numeric(outpt[i,5]))),3)
        pr<-round(1-pchisq(as.numeric(outpt[i,7]),df),4)
        pr<-ifelse(pr<0.0001,"<.0001",pr)
        outpt[i,8]<-pr
        v1<-df
        v2<-nobs-df*2
        lambda<-exp(-as.numeric(outpt[i,7])/2)
        outpt[i,9]<-round((v2/v1)*((lambda^(-2/as.numeric(nobs)))-1),3)
        prf<-round(1-pf(as.numeric(outpt[i,9]),v1,v2),4)
        prf<-ifelse(prf<0.0001,"<.0001",prf)         
        outpt[i,10]<-prf
      }
    }  
  }
  return(outpt)
}


