###############################################################################
#											            #
#                  Eggs Per Recruit Analysis                                  #
#					                                                #          
###############################################################################
#age<-c(0,1,2,3,4,5,6,7,8)               #ages
#eggs<-c(11349,23146,67846,122239,181247,236365,288478,337147,389654) #average catch weight at age
#prec<-c(0.01,0.05,0.36,0.87,0.99,1,1,1,1)                         
#pmat<-c(0,0,0.118,0.864,1,1,1,1,1)                          #partial recruitment vector
#NM<-c(1.17,0.888,0.648,0.558,0.522,0.491,0.465,0.465,0.465)                                             #natural mortality
#maxF<-5							    #maximum F for analysis
#pF<-00                                           #prop. F before spawning
#pM<-00                                           #prop. M before spawning


epr<-function(age=NULL,fecund=NULL,partial=NULL,pmat=pmat,M=NULL,pF=NULL, pM=NULL,MSP=40,plus=FALSE,oldest=NULL,maxF=2,incrF=0.0001){				
	if(is.null(age)) 
         stop ("age vector is missing") 
  	if(is.null(fecund)) 
         stop (" fecund vector is missing.") 
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
      if(any(length(age)!=c(length(age),length(fecund),length(partial),length(pmat))))
         stop("Length of vectors unequal")
      if(length(M)==1) M<-rep(M,length(age))
      data<-as.data.frame(cbind(age,fecund,partial,M,pmat,pF,pM))
	SPR<-as.data.frame(cbind(rep(NA,ceiling(maxF/incrF)+1),
                   rep(NA,ceiling(maxF/incrF)+1)))
	names(SPR)<-c("F","SPR")
      if(plus==TRUE){
                 len<-oldest-min(data$age)+1
                 if(oldest>max(data$age)){
                 pdata<-data[rep(length(data$age),times=oldest-data$age[length(data$age)]), ] 
                 pdata$age<-seq(max(data$age)+1,oldest,1)
                 data<-rbind(data,pdata)}
           }

        if(plus==FALSE) len<-max(data$age)-min(data$age)+1

#Calculate Spawning Stock Per Recruit

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
   data$SPR<-data$psb*data$SB*data$fecund*data$pmat
   SPR$SPR[i]<-sum(data$SPR)
   SPR$F[i]<-F
   F<-F+incrF
  }
  SPR$PSPR<-SPR$SPR/SPR$SPR[1]*100
   sss<-NULL
###Calculate F for percent MSP
   getF<-function(x){

   data$SB<-exp(-(data$partial*data$pF*x+data$pM*data$M))
   data$S<-cumprod(exp(-(data$partial*x+data$M)))
   data$psb[1]<-1
    for(y in 2:len)
    {
      data$psb[y]<-data$S[y-1]
    }
   data$SPR<-data$psb*data$SB*data$fecund*data$pmat
   sss<<-sum(data$SPR)
   return(((sum(data$SPR)/SPR$SPR[1]*100)-MSP)^2)
   }
  Fsp<-optimize(getF,c(0,maxF),tol=0.0000001)[1]
  ans<-NULL  
  ans<-matrix(NA,1L,2L)
  ans<-rbind(cbind(Fsp,sss)) 
  EPR<-SPR
  names(EPR)<-c("F","EPR","PEPR")
  dimnames(ans)<-list(c(paste("F at ",MSP,"%MSP",sep="")),c("F","Eggs_Per_Recruit"))  
  outpt<-list(ans,EPR);names(outpt)<-c("Reference_Point","F_vs_EPR")
  return(outpt)
}
#epr(age=age,fecund=eggs,partial=prec,
#pmat=pmat,M=NM,pF=0,pM=0,MSP=40,plus=T,maxF=4,incrF=0.001,oldest=10)


		
