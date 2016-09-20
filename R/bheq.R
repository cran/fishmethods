
bheq<-function(len=NULL,type=c(1,2),K=NULL,Linf=NULL,Lc=NULL,La=NULL,nboot=100){
  if(is.null(len)) 
         stop ("length vector does not exist.")
      if(!is.numeric(len)) 
         stop ("vector is not numeric.")
      if(is.null(K)) 
         stop ("K not specified.") 
      if(is.null(Linf)) 
         stop ("Linf not specified.") 
      if(is.null(Lc)) 
         stop ("Lc not specified.") 
      if(any(type==2) & is.null(La)) stop("La is required for method 2")
      len<-len[!is.na(len)] 
      lengths<-len[len>=Lc]
  rown<-length(type)
  results<-data.frame(method=as.character(rep("NA",rown)),
                      meanlen=as.numeric(rep(NA,rown)),
                      n=rep(as.numeric(NA),rown),
                      z=rep(as.numeric(NA),rown),
                      SE=rep(as.numeric(NA),rown),stringsAsFactors=FALSE)
  cnt=0
     if(any(type==1)){
          mean.boot1 <- function(x, i){ 
               f<-function(x){
                (K*(Linf-(sum(x)/length(x))))/((sum(x)/length(x))-Lc)
                }
               f(lengths[i])
              }
             dd<-boot(lengths,mean.boot1,R=nboot)
             mL<-mean(lengths)
          cnt<-cnt+1
         results[cnt,1]<-"BH"
         results[cnt,2:5]<-c(mL,length(lengths),dd[[1]],apply(dd$t,2,sd))
     }
  
  if(any(type==2)){
    mean.boot1 <- function(x, i){ 
      f<-function(x){
        mL<-exp(sum(log(x))/length(x))
        u<-function(z){
          d1<-((Linf-La)/(Linf-Lc))^(z/K) 
          d2<-(z*(Lc-mL)+K*(Linf-mL))/(z*(La-mL)+K*(Linf-mL))
          diff<-(d1-d2)^2
          diff
        }
        return(optimize(u,c(0.01,10),tol=0.00000001)$minimum )
      }
      f(lengths[i])
    }
    
    dd<-boot(lengths,mean.boot1,R=nboot)
    mL<-mean(lengths)
    cnt<-cnt+1
    results[cnt,1]<-"BH Bias-Corrected"
    results[cnt,2:5]<-c(mL,length(lengths),dd[[1]],apply(dd$t,2,sd))
  }
  return(results)
}


