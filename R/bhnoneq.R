bhnoneq<-function(year=NULL,mlen=NULL,ss=NULL,K=NULL,Linf=NULL,Lc=NULL,nbreaks=NULL,styrs=NULL,
    stZ=NULL, stsigma=NULL){
  if(is.null(mlen)) 
         stop ("mean length vector does not exist")
  if(is.null(year)) 
         stop ("year vector does not exist")
  if(is.null(ss)) 
         stop ("numbers vector does not exist")
      if(!is.numeric(mlen)) 
         stop ("vector is not numeric")
 if(is.null(stZ)) 
         stop ("Initial Z vector does not exist")
 if(is.null(stsigma)) 
         stop ("Initial sigma value does not exist")
      if(is.null(K)) 
         stop ("K not specified") 
      if(is.null(Linf)) 
         stop ("Linf not specified") 
      if(is.null(Lc)) 
         stop ("Lc not specified") 
      if(is.null(nbreaks)) 
         stop ("Number of mortality breaks not specified") 
     if(is.null(styrs)) 
         stop ("Starting guesses for years of mortality breaks not specified ")                  
     if(length(mlen)!=length(year))
         stop("vectors have different lengths")
  gyr<-styrs	    
  x<-as.data.frame(cbind(year,mlen,ss))
  fyr<-min(x$year);lyr<-max(x$year)
  names(x)<-c("year","mlen","m")
  nyr<-length(x[!is.na(x[,2]),2])
  count<-length(x[,1])
  mdata<-x 
  ggyr<-gyr-fyr+1
  tm<-array(0,dim=c(nbreaks,1))
  if(length(stZ)==nbreaks+1)  parms<-c(stZ,ggyr,stsigma)
  if(length(stZ)<nbreaks+1)  stop("The number of stZ values does not equal nbreak+1")
  if(length(stZ)>nbreaks+1)  stop("Too many stZ values for match to nbreak+1")

  Lpred<-NULL;results<-NULL
 Zest<-function(y){
       Z<-y[1:as.numeric(nbreaks+1)]
       sigma<-y[length(y)]
       tm[1:nbreaks,1]<-y[as.numeric(nbreaks+2):as.numeric(length(y)-1)]
       dy<-array(0,dim=c(nbreaks,count))    
       
       for(i in 1:nbreaks){
          for(j in 1:count){
              dy[i,j]<-ifelse(tm[i,1]>=j,0,j-tm[i,1])   
            }
         }
         if(nbreaks>1){
         for(i in 1:as.numeric(nbreaks-1)){
           for(j in 1:count){
             if(j>round(tm[i+1,1],0)){
                 dy[i,j]<-dy[i,j-1]
              }   
             }
           }
         }
      a<-array(0,dim=c(nbreaks+1,count))
      s<-array(0,dim=c(nbreaks+1,count))
      r<-array(0,dim=c(nbreaks+1,count))
      denom<-rep(0,count)
      numersum<-rep(0,count)
      numerator<-rep(0,count)
 
   for(m in 1:count){
       a[1,m]<-1;r[1,m]<-1  
     for(i in 1:as.numeric(nbreaks+1)){
         a[i,m]<-1;r[i,m]<-1
         if(i<as.numeric(nbreaks+1)) {s[i,m]<-1-exp(-(Z[nbreaks+2-i]+K)*dy[nbreaks+1-i,m])}
         if(i==as.numeric(nbreaks+1)){s[i,m]<-1} 
          for(j in 1:as.numeric(i-1)){
             if(i>1){
               a[i,m]<-a[i,m]*exp(-Z[nbreaks+2-j]*dy[nbreaks+1-j,m])
               r[i,m]<-r[i,m]*exp(-(Z[nbreaks+2-j]+K)*dy[nbreaks+1-j,m])
              }
            }      
         if(i<=nbreaks){denom[m]<-denom[m]+a[i,m]*((1-exp(-Z[nbreaks+2-i]*dy[nbreaks+1-i,m]))/Z[nbreaks+2-i])}
         if(i==as.numeric(nbreaks+1)) {denom[m]<-denom[m]+a[i,m]/Z[nbreaks+2-i]}
            numersum[m]<-numersum[m]+(-((1-Lc/Linf)*r[i,m]*s[i,m])/(Z[nbreaks+2-i]+K))
          }
        }
    numerator<-Linf*(denom+numersum)
    Lpred<<-numerator/denom   
    LL<--nyr*log(sigma)-sum((mdata[,3]/(2*sigma^2))*(mdata[,2]-Lpred)^2,na.rm=T)
    LL*-1
 }
   results<-optim(parms, Zest, gr = NULL,control=list(maxit=1000000,abstol=0.0000001),hessian=TRUE)
   npar<-length(parms)
   AIC<-2*results$value+2*npar
   cov<-solve(results$hessian) 
   SE<-round(sqrt(diag(cov)),3)
   tvalues<-round(results$par,2)/round(sqrt(diag(cov)),3)
   SD<-results$par[length(results$par)]
   SESD<-SE[length(results$par)]

   Zests<-results$par[1:as.numeric(nbreaks+1)]
   ZSE<-SE[1:as.numeric(nbreaks+1)]
   yrest<-results$par[as.numeric(nbreaks+2):as.numeric(length(results$par)-1)]+fyr-1
   yrSEs<-SE[as.numeric(nbreaks+2):as.numeric(length(results$par)-1)]
   vname<-c(paste("Z",seq(1,nbreaks+1,1),sep=""),paste("Y",seq(1,nbreaks,1),sep=""),"SD")
   output<-data.frame(Parameter=vname,Estimate=c(Zests,yrest,SD),SE=c(ZSE,yrSEs,SESD),t=tvalues)
   aicput<-data.frame(Parameter=c("AIC"),Estimate=c(AIC),SE=NA,t=NA)
   output<-rbind(output,aicput)
   pred<-as.data.frame(cbind(x[,1],Lpred));names(pred)<-c("year","mlen")
   output<-list(output,x,pred);names(output)<-c("results",
               "obs","pred")
   return(output)
}


