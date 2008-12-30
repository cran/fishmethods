depletM<-function(catch=NULL,effort=NULL,M=NULL,method=c("l","c"),stq=NULL,stN0=NULL){
 	if(is.null(catch)) 
         stop ("catch vector does not exist")
  	if(is.null(effort)) 
         stop ("effort vector does not exist")
  	if(length(catch)!=length(effort))
    	   stop("unequal vector lengths")
      if(any(method=="c") & is.null(stq)) stop("starting value for q is missing.")
      if(any(method=="c") & is.null(stN0)) stop("starting value for N0 is missing.")
 	x<-as.data.frame(cbind(catch,effort))
#Program
 	nsam<-length(x$catch)
       x$t<-seq(0,length(x$catch)-1,1)
 	if(length(x$catch)<3) stop("only two observations")     
  if(length(x$catch>=3)){
	#Leslie Delury model
    if(any(method=="l")){
      x$cpue<-x$catch/x$effort
      x$cpue<-ifelse(is.nan(x$cpue),0,x$cpue)
      x$S<-exp(-(x$t+0.5)*M)
     	x$CM[1]<-0
  	for(t in 2:nsam){
    	   d<-0 
     	   for(i in 1:as.numeric(t-1)){
             d<-d+x$catch[i]*exp(-(t-i)*M)
           }
        x$CM[t]<-d    
      }
 	ld<-lm(x$cpue~0+x$S+x$CM)
 	Nld<--coef(summary(ld))[1]/coef(summary(ld))[2]
 	qld<--coef(summary(ld))[2]
      predcpue<-predict(ld)
 	s2ld<-summary(ld)$sigma^2
 	SEnld<-sqrt((s2ld/qld^2)*((1/nsam)+(Nld-mean(x$CM))^2/sum((x$CM-mean(x$CM))^2)))
 	CInld<-c(Nld-qt(0.975,nsam-1)*SEnld,Nld+qt(0.975,nsam-1)*SEnld)
 	SEqld<-coef(summary(ld))[4]
 	CIqld<-c(qld-qt(0.975,nsam-1)*SEqld,qld+qt(0.975,nsam-1)*SEqld)
          lans<-NULL  
          lans$leslie.results<-matrix(NA,2L,5L)
          lans$leslie.results<-rbind(cbind(round(Nld,0),round(SEnld,1),round(CInld[1],1)
          ,round(CInld[2],1)),cbind(qld,SEqld,CIqld[1],CIqld[2]))   
          dimnames(lans$leslie.results)<-list(cbind("N0","q"),c("Estimate","SE","95% LCI","95% UCI"))
          lans$leslie.output<-as.data.frame(cbind(x$S,x$CM,predcpue,x$cpue-predcpue))
          names(lans$leslie.output)<-c("cum_M","cum_C", "Pred_CPUE","Residuals")
   }
	#catch equation model
   if(any(method=="c")){
 	x$CE[1]<-0
 	x$CE[seq(2,nsam,1)]<-cumsum(x$effort[seq(1,nsam-1,1)])
 	ce<-nls(catch~N0*exp(-(q*CE+M*(t+1)))*((q*effort)/(q*effort+M))*
         (1-exp(-(q*effort+M))),data=x,start=list(q=stq,N0=stN0))
      predC<-predict(ce) 
	Nce<-coef(summary(ce))[2]
 	SEnce<-coef(summary(ce))[4]
 	CInce<-c(Nce-qt(0.975,nsam-1)*SEnce,Nce+qt(0.975,nsam-1)*SEnce)
 	qce<-coef(summary(ce))[1]
 	SEqce<-coef(summary(ce))[3]
 	CIqce<-c(qce-qt(0.975,nsam-1)*SEqce,qce+qt(0.975,nsam-1)*SEqce)
          ceans<-NULL  
          ceans$CE.results<-matrix(NA,2L,5L)
          ceans$CE.results<-rbind(cbind(round(Nce,0),round(SEnce,1),round(CInce[1],1)
          ,round(CInce[2],1)),cbind(qce,SEqce,CIqce[1],CIqce[2]))   
          dimnames(ceans$CE.results)<-list(cbind("N0","q"),c("Estimate","SE","95% LCI","95% UCI"))
          ceans$CE.output<-as.data.frame(cbind(x$CE,predC,x$catch-predC))
          names(ceans$CE.output)<-c("cum_E","Pred_Catch","Residuals")
    }
     if(length(method)==2) return(c(lans,ceans))     
     if(length(method)==1 & any(method=="l")) return(lans)
     if(length(method)==1 & any(method=="c")) return(ceans)
 }
}
#catch=data$catch;effort=data$effort
#data<-read.csv("P:/Rwork/book/done/darter.csv")
#depletM(catch=data$catch,effort=data$effort,M=0.000,method=c("l","c"),stq=0.1,stN0=1000)


