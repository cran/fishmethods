fpc<-function(cpue1=NULL,cpue2=NULL,method=c(1,2,3,4),deletezerosets=FALSE,kapp_zeros="paired",boot_type="paired",nboot=1000,
              dint=c(1e-9,5),rint=c(1e-9,20),decimals=2,alpha=0.05){
  if(length(cpue1)!=length(cpue2)) stop("Lengths of cpue1 and cpue2 do not match.")
  if(!any(kapp_zeros %in% c("paired","ind"))) stop("unknown argument for kapp_zeros")
  if(!any(deletezerosets %in% c(TRUE,FALSE))) stop("unknown argument for deletezerosets")
  if(!any(boot_type %in% c("paired","unpaired"))) stop("unknown argument for boot_type")
  combdata<-data.frame(cpue1=cpue1,cpue2=cpue2)
  #remove 0,0 rows 
  combdata<-combdata[!apply(combdata, 1, function(x) all(x<=0)),]
  if(deletezerosets==TRUE) combdata<-combdata[!apply(combdata, 1, function(x) any(x<=0) ),] 
  randomRows <- function(file,n){
         return(file[sample(nrow(file),n,replace=TRUE),])
       }
  npairs<-length(combdata[,1])
  z<-abs(qnorm(alpha/2))
 ################## Ratio of Means ############################
  
 if(any(method==1)){
  M1data<-combdata
  M1<-mean(M1data$cpue1)/mean(M1data$cpue2)
  mcpue2<-mean(M1data$cpue2)
  #from Thompson 2002 page 70
  vr<-sum((M1data$cpue1-M1*M1data$cpue2)^2)/(npairs-1)
  VARr<-(1/mcpue2^2)*(vr/npairs)
  SEM1U<-sqrt(VARr)
  U_LCI_M1<-M1-z*SEM1U
  U_UCI_M1<-M1+z*SEM1U
  #Variance Jackknife
  if(boot_type=="paired"){ 
    tempM1<-NULL
    for(jack in 1:length(M1data[,1])){
      jackdata<-M1data[-jack,]  
      tempM1[jack]<-mean(jackdata$cpue1)/mean(jackdata$cpue2)
     }  
    pseudo<-npairs*M2-(npairs-1)*tempM1
    SEM1J<-sd(pseudo)/sqrt(npairs)    
    J_LCI_M1<-M1-z*SEM1J
    J_UCI_M1<-M1+z*SEM1J
    outsM1<-jackknife(1:length(M1data[,1]),function(x,M1data){(sum(M1data[x,1])/(npairs-1))/(sum(M1data[x,2])/(npairs-1))},M1data)
    
  }
  if(boot_type=="unpaired"){
    SEM1J<-NA 
    J_LCI_M1<-NA
    J_UCI_M1<-NA
  }
  #Bootstrap var
  tempM1<-NULL
  for(boots in 1:nboot){
    if(boot_type=="unpaired"){
      y<-sample(M1data[,1],length(M1data[,1]),replace=TRUE)
      x<-sample(M1data[,2],length(M1data[,2]),replace=TRUE)
    }
    if(boot_type=="paired"){
      randomsample<-randomRows(M1data,length(M1data[,1]))
      y<-randomsample[,1]
      x<-randomsample[,2]
    }
    tempM1[boots]<-(sum(y)/npairs)/(sum(x)/npairs)
  }
  SEM1B<-sd(tempM1) 
  B_CIM1<-quantile(tempM1,probs=c(alpha/2,1-(alpha/2)))
 }
################# Randomized Block ANOVA ##############################
  if(any(method==2)){
    M2data<-data.frame(group=c(rep(1,npairs),rep(2,npairs)),haul=c(1:npairs,1:npairs),
                      cpue=c(combdata$cpue1,combdata$cpue2),logcpue=c(log(combdata$cpue1+1),log(combdata$cpue2+1)))
    options(contrasts=c("contr.sum","contr.poly"))
    outsM2<-summary(glm(logcpue~factor(group)+factor(haul),data=M2data))$coef[1:2,1:2]
    M2<-exp(2*outsM2[2,1]*(1+0.5*outsM2[2,2]^2))
    SEM2U<-outsM2[2,2]
    U_LCI_M2<-exp(2*outsM2[2,1]-2*z*SEM2U)
    U_UCI_M2<-exp(2*outsM2[2,1]+2*z*SEM2U)

    #Bootstrap
    tempM2<-NULL
    for(boots in 1:nboot){
      if(boot_type=="unpaired"){
        y<-sample(combdata[,1],length(combdata[,1]),replace=TRUE)
        x<-sample(combdata[,2],length(combdata[,2]),replace=TRUE)
      }
      if(boot_type=="paired"){
        randomsample<-randomRows(combdata,length(combdata[,1]))
        y<-randomsample[,1]
        x<-randomsample[,2]
      }
      M2data<-data.frame(group=c(rep(1,npairs),rep(2,npairs)),haul=c(1:npairs,1:npairs),
                         cpue=c(y,x),
                         logcpue=c(log(y+1),log(x+1)))
      outsM2<-summary(glm(logcpue~factor(group)+factor(haul),data=M2data))$coef[1:2,1:2]
      tempM2[boots]<-exp(2*outsM2[2,1]*(1+0.5*outsM2[2,2]^2))
    }  
    SEM2B<-sd(tempM2)    
    B_CIM2<-quantile(tempM2,probs=c(alpha/2,1-(alpha/2)))  
    #jacknife
   if(boot_type=="paired"){ 
    tempM2<-NULL
    for(jack in 1:length(combdata[,1])){
      jackdata<-combdata[-jack,]  
      M2data<-data.frame(group=c(rep(1,npairs-1),rep(2,npairs-1)),haul=c(1:c(npairs-1),1:c(npairs-1)),
                         cpue=c(jackdata$cpue1,jackdata$cpue2),
                         logcpue=c(log(jackdata$cpue1+1),log(jackdata$cpue2+1)))
      outsM2<-summary(glm(logcpue~factor(group)+factor(haul),data=M2data))$coef[1:2,1:2]
      tempM2[jack]<-exp(2*outsM2[2,1]*(1+0.5*outsM2[2,2]^2))
    }  
    pseudo<-npairs*M2-(npairs-1)*tempM2
    SEM2J<-sd(pseudo)/sqrt(npairs)    
    J_LCI_M2<-M2-z*SEM2J
    J_UCI_M2<-M2+z*SEM2J
   }
   if(boot_type=="unpaired"){ 
    SEM2J<-NA    
    J_LCI_M2<-NA
    J_UCI_M2<-NA
   }
  }#Method 2
  ############### Method 3 ################################
  if(any(method==3)){
    M3data<-data.frame(group=c(rep(1,npairs),rep(2,npairs)),haul=c(1:npairs,1:npairs),
                       cpue=c(combdata$cpue1,combdata$cpue2),logcpue=c(log(combdata$cpue1+1),log(combdata$cpue2+1)))
   options(contrasts=c("contr.sum","contr.poly"))
   outsM3<-summary(glm(logcpue~factor(group),data=M3data))$coef[1:2,1:2]
   M3<-exp(outsM3[1,1]+outsM3[2,1]*(1+0.5*outsM3[2,2]^2))/exp(outsM3[1,1]-outsM3[2,1]*(1+0.5*outsM3[2,2]^2))
   SEM3<-outsM3[2,2]
   U_UCI_M3<-exp(2*outsM3[2,1]+2*z*SEM3)
   U_LCI_M3<-exp(2*outsM3[2,1]-2*z*SEM3)
   #Bootstrap
   tempM3<-NULL
   for(boots in 1:nboot){
     if(boot_type=="unpaired"){
       y<-sample(combdata[,1],length(combdata[,1]),replace=TRUE)
       x<-sample(combdata[,2],length(combdata[,2]),replace=TRUE)
     }
     if(boot_type=="paired"){
       randomsample<-randomRows(combdata,length(combdata[,1]))
       y<-randomsample[,1]
       x<-randomsample[,2]
     }
     M3data<-data.frame(group=c(rep(1,npairs),rep(2,npairs)),haul=c(1:npairs,1:npairs),
                        cpue=c(y,x),
                        logcpue=c(log(y+1),log(x+1)))
     outsM3<-summary(glm(logcpue~factor(group),data=M3data))$coef[1:2,1:2]
     tempM3[boots]<-exp(2*outsM3[2,1]*(1+0.5*outsM3[2,2]^2))
   }  
   SEM3B<-sd(tempM3)    
   B_CIM3<-quantile(tempM3,probs=c(alpha/2,1-(alpha/2)))  
   #jacknife
   if(boot_type=="paired"){
   tempM3<-NULL
   for(jack in 1:length(combdata[,1])){
     jackdata<-combdata[-jack,]  
     M3data<-data.frame(group=c(rep(1,npairs-1),rep(2,npairs-1)),haul=c(1:c(npairs-1),1:c(npairs-1)),
                        cpue=c(jackdata$cpue1,jackdata$cpue2),
                        logcpue=c(log(jackdata$cpue1+1),log(jackdata$cpue2+1)))
     outsM3<-summary(glm(logcpue~factor(group),data=M3data))$coef[1:2,1:2]
     tempM3[jack]<-exp(2*outsM3[2,1]*(1+0.5*outsM3[2,2]^2))
   }  
   pseudo<-npairs*M3-(npairs-1)*tempM3
   SEM3J<-sd(pseudo)/sqrt(npairs)    
   J_LCI_M3<-M3-z*SEM3J
   J_UCI_M3<-M3+z*SEM3J
   }
   if(boot_type=="unpaired"){
     SEM3J<-NA    
     J_LCI_M3<-NA
     J_UCI_M3<-NA
   } 
  }
 
  ############### Method 4 #################################
  if(any(method==4)){
    if(kapp_zeros=="ind"){
      yold<-combdata$cpue1
      yold<-yold[yold>0]
      xold<-combdata$cpue2
      xold<-xold[xold>0]
      omeancpue1<-mean(yold)
      omeancpue2<-mean(xold)
      orign<-length(xold)
      origm<-length(yold)
    }
    if(kapp_zeros=="paired"){
      tempM4<-combdata[!apply(combdata, 1, function(x) any(x<=0) ),]
      omeancpue1<-mean(tempM4[,1])
      omeancpue2<-mean(tempM4[,2])
      yold<-tempM4[,1]
      xold<-tempM4[,2]
      orign<-length(xold)
      origm<-length(yold)
    }
    # Original Estimate
    getd<-function(d,a){
      n<-length(x)
      m<-length(y)
      xd<-x^d
      yd<-y^d
      x2d<-x^(2*d)
      y2d<-y^(2*d)
      t<-(m/((n+m)^2))*(sum(xd)^2)-(m/(n+m))*sum(x2d)
      s<-((m-n)/((n+m)^2))*sum(xd)*sum(yd)
      p<-(n/(n+m))*sum(y2d)-(n/((n+m)^2))*sum(yd)^2
      q<-(-s+sqrt((s^2)-4*p*t))/2*p
      v<-(sum(xd)+q*sum(yd))/(n+m)
      w<-(1/(n+m))*(sum((xd-v)^2)+sum(((q*yd)-v)^2))
      u<-sum(log(x))+sum(log(y))
      value<-((n+m)/d)+u-(1/w)*(sum(x2d*log(x))-v*sum(log(x)*xd)+
                                  (q^2)*sum(y2d*log(y)))-v*q*sum(yd*log(y))
      return(value-a)
    }
    solver<-function(r1){
      pred_xdj<-NULL;pred_ydj<-NULL
      n<-length(x)
      m<-length(y)
      for(jj in 1:length(x)) pred_xdj[jj]<-1/(n+m-1)*(sum(x[-jj]^d)+sum((r1^d)*(y^d)))
      for(jj in 1:length(y)) pred_ydj[jj]<-1/(n+m-1)*(sum(x^d)+sum((r1^d)*(y[-jj]^d)))
      Sa<-sum(((x^d)-pred_xdj)^2)+sum((((r1^d)*(y^d))-pred_ydj)^2)
      pred_xdj2<-NULL;pred_ydj2<-NULL
      for(jj in 1:length(x)) pred_xdj2[jj]<-1/(n-1)*sum(x[-jj]^d)
      for(jj in 1:length(y)) pred_ydj2[jj]<-1/(m-1)*sum((r1^d)*(y[-jj]^d))
      Sb<-sum(((x^d)-pred_xdj2)^2)+sum((((r1^d)*(y^d))-pred_ydj2)^2)
      return(Sa-Sb)
    }
    y<-yold
    x<-xold 
    d<-uniroot(getd, c(dint[1], dint[2]), extendInt="yes",tol = 1e-9, a = 0)[[1]]
    M4<-1/optimize(solver,c(rint[1],rint[2]))$minimum
    
    #Jackknife
    SEM4J<-NA
    J_LCI_M4<-NA
    J_UCI_M4<-NA
    if(kapp_zeros == "paired" & boot_type=="paired"){
      jackdata<-NULL
      for(jack in 1:length(y)){
        y<-tempM4[-jack,1]
        x<-tempM4[-jack,2]
        d<-uniroot(getd, c(dint[1], dint[2]), extendInt="yes",tol = 1e-9, a = 0)[[1]]
        jackdata[jack]<-1/optimize(solver,c(rint[1],rint[2]))$minimum
      }
      pseudo<-length(tempM4[,1])*M4-(length(tempM4[,1])-1)*jackdata
      SEM4J<-sd(pseudo)/sqrt(length(tempM4[,1]))    
      J_LCI_M4<-M4-z*SEM4J
      J_UCI_M4<-M4+z*SEM4J
    }  
    
    #Bootstrap
    M4store<-NULL
    for(boots in 1:nboot){
     if(kapp_zeros %in% c("paired","ind") & boot_type=="unpaired"){
       y<-sample(yold,length(yold),replace=TRUE)
       x<-sample(xold,length(xold),replace=TRUE)
     }
     if(kapp_zeros == "paired" & boot_type=="paired"){
       randomsample<-randomRows(tempM4,length(tempM4[,1]))
       y<-randomsample[,1]
       x<-randomsample[,2]
     }
     d<-uniroot(getd, c(dint[1], dint[2]), extendInt="yes",tol = 1e-9, a = 0)[[1]]
     M4store[boots]<-1/optimize(solver,c(rint[1],rint[2]))$minimum
   }#bootstrap
    SEM4B<-sd(M4store)
    B_CI4<-quantile(M4store,probs=c(alpha/2,1-(alpha/2)))
    if(kapp_zeros=="ind" & boot_type=="paired"){
      SEM4B<-NA
      B_CI4<-c(NA,NA)
    }
  }# method 4

  mrow<-length(method)
  outtable<-matrix(nrow=mrow,ncol=15)
  labs<-round((1-alpha)*100,0)
  colnames(outtable)<-c("method","n1","n2","mean cpue1","mean cpue2","FPC","U_SE","Jack_SE","Boot_SE",
                        paste("U_",labs,"%_LCI",sep=""),paste("U_",labs,"%_UCI",sep=""),
                        paste("Jack_",labs,"%_LCI",sep=""),paste("Jack_",labs,"%_UCI",sep=""),
                        paste("Boot_",labs,"%_LCI",sep=""),paste("Boot_",labs,"%_UCI",sep=""))
  outtable<-as.data.frame(outtable)
  cnt<-0
  if(any(method==1)){
    cnt<-cnt+1
    outtable[cnt,1]<-"Ratio of Means"
    outtable[cnt,2]<-npairs
    outtable[cnt,3]<-npairs
    outtable[cnt,4]<-round(mean(M1data[,1]),decimals)
    outtable[cnt,5]<-round(mean(M1data[,2]),decimals)
    outtable[cnt,6]<-round(M1,decimals)
    outtable[cnt,7]<-round(SEM1U,decimals)
    outtable[cnt,8]<-round(SEM1J,decimals)
    outtable[cnt,9]<-round(SEM1B,decimals)
    outtable[cnt,10]<-round(U_LCI_M1,decimals)
    outtable[cnt,11]<-round(U_UCI_M1,decimals)
    outtable[cnt,12]<-round(J_LCI_M1,decimals)
    outtable[cnt,13]<-round(J_LCI_M1,decimals)
    outtable[cnt,14]<-round(as.numeric(B_CIM1[1]),decimals)
    outtable[cnt,15]<-round(as.numeric(B_CIM1[2]),decimals)
  }
 if(any(method==2)){
   cnt<-cnt+1
   outtable[cnt,1]<-"Randomized Block ANOVA"
   outtable[cnt,2]<-npairs
   outtable[cnt,3]<-npairs
   outtable[cnt,4]<-round(mean(combdata[,1]),decimals)
   outtable[cnt,5]<-round(mean(combdata[,2]),decimals)
   outtable[cnt,6]<-round(M2,decimals)
   outtable[cnt,7]<-NA
   outtable[cnt,8]<-round(SEM2J,decimals)
   outtable[cnt,9]<-round(SEM2B,decimals)
   outtable[cnt,10]<-round(U_LCI_M2,decimals)
   outtable[cnt,11]<-round(U_UCI_M2,decimals)
   outtable[cnt,12]<-round(J_LCI_M2,decimals)
   outtable[cnt,13]<-round(J_UCI_M2,decimals)
   outtable[cnt,14]<-round(as.numeric(B_CIM2[1]),decimals)
   outtable[cnt,15]<-round(as.numeric(B_CIM2[2]),decimals)
 }
 if(any(method==3)){
   cnt<-cnt+1
   outtable[cnt,1]<-"Multiplicative Model"
   outtable[cnt,2]<-npairs
   outtable[cnt,3]<-npairs
   outtable[cnt,4]<-round(mean(combdata[,1]),decimals)
   outtable[cnt,5]<-round(mean(combdata[,2]),decimals)
   outtable[cnt,6]<-round(M3,decimals)
   outtable[cnt,7]<-NA
   outtable[cnt,8]<-round(SEM3J,decimals)
   outtable[cnt,9]<-round(SEM3B,decimals)
   outtable[cnt,10]<-round(U_LCI_M3,decimals)
   outtable[cnt,11]<-round(U_UCI_M3,decimals)
   outtable[cnt,12]<-round(J_LCI_M3,decimals)
   outtable[cnt,13]<-round(J_UCI_M3,decimals)
   outtable[cnt,14]<-round(as.numeric(B_CIM3[1]),decimals)
   outtable[cnt,15]<-round(as.numeric(B_CIM3[2]),decimals)
 }
 if(any(method==4)){
   cnt<-cnt+1
   outtable[cnt,1]<-"Kappenman"
   outtable[cnt,2]<-origm
   outtable[cnt,3]<-orign
   outtable[cnt,4]<-round(omeancpue1,decimals)
   outtable[cnt,5]<-round(omeancpue2,decimals)
   outtable[cnt,6]<-round(M4,decimals)
   outtable[cnt,7]<-NA
   outtable[cnt,8]<-round(SEM4J,decimals)
   outtable[cnt,9]<-round(SEM4B,decimals)
   outtable[cnt,10]<-NA
   outtable[cnt,11]<-NA
   outtable[cnt,12]<-round(J_LCI_M4,decimals)
   outtable[cnt,13]<-round(J_UCI_M4,decimals)
   outtable[cnt,14]<-round(as.numeric(B_CI4[1]),decimals)
   outtable[cnt,15]<-round(as.numeric(B_CI4[2]),decimals)
 }
 return(outtable)
}
   