sbpr<-function(age=NULL,ssbwgt=NULL,partial=NULL,pmat=NULL,M=NULL,pF=NULL,pM=NULL,
               plus=FALSE,oldest=NULL,maxF=2,incrF=0.0001,options=c(1,2,3,4),MSP=NULL,SSBPR=NULL,
               Fsol=NULL,graph=TRUE){		
  # 8/21/2023 added options arguments and procedures
	if(is.null(age)) 
         stop ("age vector is missing") 
  	if(is.null(ssbwgt)) 
         stop (" ssbwgt vector is missing.") 
     	if(is.null(partial)) 
         stop ("partial recruitment vector is missing.")
     	if(is.null(pmat)) 
         stop ("pmat vector is missing.")
  	if(is.null(M)) 
         stop ("M value or vector is missing")
      if(is.null(pF))
           stop ("pF value is missing.")
      if(is.null(pM))
           stop ("pM value is missing.")
      if(plus==TRUE & is.null(oldest)) stop("oldest must be specified for plus group calculation.")    
      if(any(length(age)!=c(length(age),length(ssbwgt),length(partial),length(pmat))))
         stop("Length of vectors unequal")
      if(length(M)==1) M<-rep(M,length(age))
      data<-as.data.frame(cbind(age,ssbwgt,partial,M,pmat,pF,pM))
      if(plus==TRUE){
                 len<-oldest-min(data$age)+1
                 if(oldest>max(data$age)){
                 pdata<-data[rep(length(data$age),times=oldest-data$age[length(data$age)]),] 
                 pdata$age<-seq(max(data$age)+1,oldest,1)
                 data<-rbind(data,pdata)}
           }
      if(plus==FALSE) len<-max(data$age)-min(data$age)+1
      if(length(which(!options %in% c(1,2,3,4)))>0) stop("invalid option number") 
      
#Calculate Spawning Stock Per Recruit Increment Matrix
if(any(options==1)){
  SPR<-as.data.frame(cbind(rep(NA,ceiling(maxF/incrF)+1),
                           rep(NA,ceiling(maxF/incrF)+1)))
  names(SPR)<-c("F","SPR")
F<-0
for (i in 1:length(SPR$F))
  {
   data$SB<-exp(-(data$partial*data$pF*F+data$pM*data$M))
   data$S<-cumprod(exp(-(data$partial*F+data$M)))
   data$psb[1]<-1
  for(y in 2:len)
   {
    data$psb[y]<-data$S[y-1]
   }
   data$SPR<-data$psb*data$SB*data$ssbwgt*data$pmat
   SPR[i,2]<-sum(data$SPR)
   SPR[i,1]<-F
   F<-F+incrF
  }
  SPR[,3]<-SPR[,2]/SPR[1,2]*100
  names(SPR)<-c("F","SSBPR","PSSBPR")
}
#Find SBPR for single F
if(any(options==2)){
  if(is.null(Fsol)) stop ("Fsol value is missing for option 2.")
	    data$SB<-exp(-(data$partial*data$pF*Fsol+data$pM*data$M))
	    data$S<-cumprod(exp(-(data$partial*Fsol+data$M)))
	    data$psb[1]<-1
	    for(y in 2:len)
	    {
	      data$psb[y]<-data$S[y-1]
	    }
	    data$SPR<-data$psb*data$SB*data$ssbwgt*data$pmat
	    SPR2<-data.frame(Fsol=Fsol,SSBPR=sum(data$SPR))
} 
	

###Calculate F for percent MSP
if(any(options==3)){
  if(is.null(MSP)) stop ("MSP value is missing for option 3.")
  data$SB<-exp(-(data$partial*data$pF*0+data$pM*data$M))
  data$S<-cumprod(exp(-(data$partial*0+data$M)))
  data$psb[1]<-1
  for(y in 2:len)
  {
    data$psb[y]<-data$S[y-1]
  }
  data$SPR<-data$psb*data$SB*data$ssbwgt*data$pmat
  SPRF0<-sum(data$SPR)
  
   getF<-function(x){
   data$SB<-exp(-(data$partial*data$pF*x+data$pM*data$M))
   data$S<-cumprod(exp(-(data$partial*x+data$M)))
   data$psb[1]<-1
    for(y in 2:len)
    {
      data$psb[y]<-data$S[y-1]
    }
   data$SPR<-data$psb*data$SB*data$ssbwgt*data$pmat
   return(((sum(data$SPR)/SPRF0*100)-MSP)^2)
   }
  Fmsp<-optimize(getF,c(0,maxF),tol=0.00000001)[1]
  SPR3<-data.frame(F=as.numeric(Fmsp),PSSBPR=MSP)
}
	###Find F for SBR
	if(any(options==4)){
	  if(is.null(SSBPR)) stop ("SSBPR value is missing for option 4.")

	  getF<-function(x){
	    data$SB<-exp(-(data$partial*data$pF*x+data$pM*data$M))
	    data$S<-cumprod(exp(-(data$partial*x+data$M)))
	    data$psb[1]<-1
	    for(y in 2:len)
	    {
	      data$psb[y]<-data$S[y-1]
	    }
	    data$SPR<-data$psb*data$SB*data$ssbwgt*data$pmat
	    return((sum(data$SPR)-SSBPR)^2)
	  }
	  Fsbr<-optimize(getF,c(0,maxF),tol=0.00000001)[1]
	  SPR4<-data.frame(F=as.numeric(Fsbr),SSBPR=SSBPR)
	}
    ans<-NULL
    labels<-NULL
  if(any(options==2)){
    ans<-c(ans,list(SPR2)) 
    labels<-c(labels,"SSBPR_at_Fsol")
    }
  if(any(options==3)){
    ans<-c(ans,list(SPR3))    
    labels<-c(labels,"F_at_MSP")
  }
  if(any(options==4)){
    ans<-c(ans,list(SPR4))    
    labels<-c(labels,"F_at_SSBPR")
  }
  if(any(options==1)){
    ans<-c(ans,list(SPR))
    labels<-c(labels,"F_vs_SSBPR")
    }
  names(ans)<-labels  
  if(any(options==1) & graph==TRUE){
    par(mfrow=c(1,2))
    plot(SPR[,2]~SPR[,1],ylab="SSBPR",xlab="F",type="l")
    plot(SPR[,3]~SPR[,1],ylab="% Max SSBPR",xlab="F",type="l")
  }
  return(ans)
}