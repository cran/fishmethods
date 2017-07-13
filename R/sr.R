sr<-function(recruits=NULL,stock=NULL,model=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),select=1,
    initial=list(RA=NULL,RB=NULL,Rrho=NULL,
                 BHA=NULL,BHB=NULL,BHrho=NULL,
                 SHA=NULL,SHB=NULL,SHC=NULL, 
                 DSA=NULL,DSB=NULL,DSC=NULL,
                 MYA=NULL,MYB=NULL,MYC=NULL),
    control=list(maxit=10000),plot=FALSE){
if(is.null(recruits)) stop("recruits vector is empty")
if(is.null(stock)) stop("stock vector is empty")
if(length(recruits)!=length(stock)) stop("recruits and stock vectors are not the same length")
datar<-as.data.frame(cbind(recruits=recruits,stock=stock))
datar<-datar[!is.na(datar[,1]),]
datar<-datar[!is.na(datar[,2]),]
datar<-datar[order(datar[,2]),]
n<-length(datar$recruits)
################## Model 0 ##################################################
if(any(model==0)){ #Density-Independent
 parms<-1
 di<-function(parms){
     vpred<-parms[1]*datar$stock        
     sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
     NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+
          1/(2*sigma2)*sum((datar$recruits-vpred)^2)
     NLL
   }
  model0<-nlm(di,parms, hessian = TRUE)
    if(class(model0)!="try-error"){
     K<-2
     model0$pred<-vpred<-model0$estimate*datar$stock
     model0$sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
     model0$AICC<-2*model0$minimum+2*K*(n/(n-K-1))
     model0$resids<-datar$recruits-vpred
     model0$convergence<-0
   }
 }  
##################  Model 1 #################################################
if(any(model==1)){ #Ricker normal errors
   if(select==1){
    temp<-lm(log(datar$recruits/datar$stock)~datar$stock)
    parms<-c(as.numeric(exp(coef(temp)[1])),as.numeric(-coef(temp)[2]))
   }
   if(select==2){
     if(is.null(initial$RA)) stop("starting value for RA is missing")
     if(is.null(initial$RB)) stop("starting value for RB is missing")
        parms<-c(initial$RA,initial$RB)
   }
    ricknorm<-function(parms){
		  A<-parms[1]
      B<-parms[2]
      vpred<-A*datar$stock*exp(-B*datar$stock)  
      sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
      NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+
               1/(2*sigma2)*sum((datar$recruits-vpred)^2)
      NLL
     }  
    model1 <- try(optim(parms,ricknorm,control=control),silent=FALSE)
 if(class(model1)!="try-error"){
     K<-length(model1$par)+1 #added 1 for sigma2
     model1$hessian<-hessian(ricknorm,model1$par)      
     model1$pred<-vpred<-model1$par[1]*datar$stock*exp(-model1$par[2]*datar$stock)
     model1$sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
     model1$AICC<-2*model1$value+2*K*(n/(n-K-1))
     model1$resids<-datar$recruits-vpred   
  }
}
################ Model 2 ##############################################
if(any(model==2)){ #Ricker log-normal errors
 if(select==1){
   temp<-lm(log(datar$recruits/datar$stock)~datar$stock)
    parms<-c(as.numeric(exp(coef(temp)[1])),as.numeric(-coef(temp)[2]))
   
     }
   if(select==2){
     if(is.null(initial$RA)) stop("starting value for parameter RA is missing")
     if(is.null(initial$RB)) stop("starting value for parameter RB is missing")
     parms<-c(initial$RA,initial$RB)
   }
   ricklog<-function(parms){
		    A<-parms[1]
        B<-parms[2]
        vpred<-A*datar$stock*exp(-B*datar$stock)        
		    vpred<-ifelse(vpred<=0,0.0001,vpred)
		    sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
        NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+sum(log(datar$recruits))+
         1/(2*sigma2)*sum((log(datar$recruits)-log(vpred)+sigma2/2)^2)
        NLL      
     }  
    model2 <- try(optim(parms,ricklog,control=control),silent=FALSE)
if(class(model2)!="try-error"){
     K<-length(model2$par)+1 #add 1 for sigma2  
     model2$hessian<-hessian(ricklog,model2$par)
     model2$pred<-vpred<-model2$par[1]*datar$stock*exp(-model2$par[2]*datar$stock)
     model2$sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
     model2$AICC<-2*model2$value+2*K*(n/(n-K-1))
	   model2$resids<-log(datar$recruits)-log(vpred)
    }
  } # Ricker log-normal
############################ Model 3 ##################################
if(any(model==3)){ # Ricker correlated normal errors
if(select==1){
    temp<-lm(log(datar$recruits/datar$stock)~datar$stock)
    parms<-c(as.numeric(exp(coef(temp)[1])),as.numeric(-coef(temp)[2]))
    tres<-datar$recruits-(parms[1]*datar$stock*exp(-parms[2]*datar$stock))
    ar1<-arima(tres,order=c(1,0,0),include.mean=FALSE)
    parms<-c(parms,as.numeric(ar1$coef))
 }
 if(select==2){
     if(is.null(initial$RA)) stop("starting value for parameter RA is missing")
     if(is.null(initial$RB)) stop("starting value for parameter RB is missing")
     if(is.null(initial$Rrho)) stop("starting value for parameter Rrho is missing")
      parms<-c(initial$RA,initial$RB,initial$Rrho)
  }
    ricknormcorr<-function(parms){
		  A<-parms[1]
      B<-parms[2]              
		  phi<-parms[3]
      vpred<-A*datar$stock*exp(-B*datar$stock)  
      phi<-ifelse(phi>1,1,ifelse(phi<=-1,-1,phi))
      res<-datar$recruits-vpred
      es<-c(res[1:c(length(res)-1)]*phi)
      sigma2w<-sum((res[-1]-es)^2)/length(res[-1])
        sump1<-0
      for(k in 2:n) sump1<-sump1+(datar$recruits[k]-phi*datar$recruits[k-1]-vpred[k]+phi*vpred[k-1])^2
       NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2w))-0.5*log(1-phi^2)+
       1/(2*sigma2w)*sump1+((1-phi^2)/(2*sigma2w))*(datar$recruits[1]-vpred[1])^2           
       NLL  
    }  
    model3 <- try(optim(parms,ricknormcorr,control=control),silent=TRUE)
  if(class(model3)!="try-error"){
      K<-length(model3$par)+1
      model3$hessian<-hessian(ricknormcorr,model3$par)
      model3$pred<-vpred<-model3$par[1]*datar$stock*exp(-model3$par[2]*datar$stock)
      res<-datar$recruits-vpred
      es<-c(res[1:c(length(res)-1)]*model3$par[3])
      sigma2w<-sum((res[-1]-es)^2)/length(res[-1])
      model3$sigma2w<-sigma2w
      model3$AICC<-2*model3$value+2*K*(n/(n-K-1))
	    model3$resids<-datar$recruits-vpred
   }
}#Ricker normal correlated
############################### Model 4 #################################
if(any(model==4)){ # Ricker correlated log-normal errors
if(select==1){
     temp<-lm(log(datar$recruits/datar$stock)~datar$stock)
     parms<-c(as.numeric(exp(coef(temp)[1])),as.numeric(-coef(temp)[2]))
     tres<-log(datar$recruits)-log(parms[1]*datar$stock*exp(-parms[2]*datar$stock))
     ar1<-arima(tres,order=c(1,0,0),include.mean=FALSE)
     parms<-c(parms,as.numeric(ar1$coef))
  }
   if(select==2){
    if(is.null(initial$RA)) stop("starting value for parameter RA is missing")
    if(is.null(initial$RB)) stop("starting value for parameter RB is missing")
    if(is.null(initial$Rrho)) stop("starting value for parameter Rrho is missing")
     parms<-c(initial$RA,initial$RB,initial$Rrho)
   }
    ricklogcorr<-function(parms){
		   A<-parms[1]
       B<-parms[2]
       phi<-parms[3] 
       phi<-ifelse(phi>1,1,ifelse(phi<=-1,-1,phi))
       vpred<-A*datar$stock*exp(-B*datar$stock)     
		   vpred<-ifelse(vpred<=0,0.0001,vpred)
		   res<-log(datar$recruits)-log(vpred)
		   es<-c(res[1:c(length(res)-1)]*phi)
		   sigma2w<-sum((res[-1]-es)^2)/length(res[-1])  
       sump1<-0
       for(k in 2:n) sump1<-sump1+(log(datar$recruits[k])-phi*log(datar$recruits[k-1])-log(vpred[k])+phi*log(vpred[k-1])+(1-phi)*sigma2w/2)^2
         NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2w))+sum(log(datar$recruits))-0.5*log(1-phi^2)+
         1/(2*sigma2w)*sump1+((1-phi^2)/(2*sigma2w))*(log(datar$recruits[1])-log(vpred[1])+sigma2w/2)^2           
        NLL  
    }  
    model4 <- try(optim(parms,ricklogcorr,control=control),silent=TRUE)
 if(class(model4)!="try-error"){
      K<-length(model4$par)+1
      model4$hessian<-hessian(ricklogcorr,model4$par)
      model4$pred<-vpred<-model4$par[1]*datar$stock*exp(-model4$par[2]*datar$stock)
      res<-log(datar$recruits)-log(vpred)
      es<-c(res[1:c(length(res)-1)]*model4$par[3])
      model4$sigma2w<-sum((res[-1]-es)^2)/length(res[-1])  
      model4$AICC<-2*model4$value+2*K*(n/(n-K-1))
      model4$resids<-log(datar$recruits)-log(vpred)
    }
  
}#Ricker log-normal correlated

##################  Model 5 #################################################
if(any(model==5)){ #Beverton-Holt normal errors
   if(select==1){
     index<-which(datar$recruits>0)[1:2]
     parms<-c(sum(datar$recruits[index])/sum(datar$stock[index]),mean(datar$recruits[c(length(datar$recruits)/2):length(datar$recruits)]))
   }
   if(select==2){
     if(is.null(initial$BHA)) stop("starting value for parameter BHA is missing")
     if(is.null(initial$BHB)) stop("starting value for parameter BHB is missing")
          parms<-c(initial$BHA,initial$BHB)
   }
    bhnorm<-function(parms){
		  A<-parms[1]
      B<-parms[2]
      vpred<-(A*datar$stock)/(1+(A*datar$stock)/B)  
      sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
		  NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+
         1/(2*sigma2)*sum((datar$recruits-vpred)^2)
      NLL
    }  
    model5 <- try(optim(parms,bhnorm,control=control),silent=FALSE)
 if(class(model5)!="try-error"){
    K<-length(model5$par)+1 # add 1 for sigma2
    model5$hessian<-hessian(bhnorm,model5$par)
    model5$pred<-vpred<-(model5$par[1]*datar$stock)/(1+(model5$par[1]*datar$stock)/model5$par[2])
    model5$sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    model5$AICC<-2*model5$value+2*K*(n/(n-K-1))
    model5$resids<-datar$recruits-vpred
   }
}
################ Model 6 ##############################################
if(any(model==6)){ #Beverton-Holt log-normal errors
 if(select==1){
   index<-which(datar$recruits>0)[1:2]
    parms<-c(sum(datar$recruits[index])/sum(datar$stock[index]),mean(datar$recruits[c(length(datar$recruits)/2):length(datar$recruits)]))
   }
   if(select==2){
     if(is.null(initial$BHA)) stop("starting value for parameter BHA is missing")
     if(is.null(initial$BHB)) stop("starting value for parameter BHB is missing")
       parms<-c(initial$BHA,initial$BHB)
   }
    bhlog<-function(parms){
		  A<-parms[1]
      B<-parms[2]
      vpred<-(A*datar$stock)/(1+(A*datar$stock)/B)
      vpred<-ifelse(vpred<=0,0.0001,vpred)
      sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
      NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+sum(log(datar$recruits))+
           1/(2*sigma2)*sum((log(datar$recruits)-log(vpred)+sigma2/2)^2)
      NLL      
    }  
    model6 <- try(optim(parms,bhlog,control=control),silent=FALSE)
    if(class(model6)!="try-error"){
      K<-length(model6$par)+1 #add 1 for sigma
      model6$hessian<-hessian(bhlog,model6$par)    
      model6$pred<-vpred<-(model6$par[1]*datar$stock)/(1+(model6$par[1]*datar$stock)/model6$par[2])
      model6$sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
      model6$AICC<-2*model6$value+2*K*(n/(n-K-1))
      model6$resids<-log(datar$recruits)-log(vpred)
    }
  } # BHlog-normal
############################ Model 7 ##################################
if(any(model==7)){ # Beverton-Holt correlated normal errors
if(select==1){
  index<-which(datar$recruits>0)[1:2]
     parms<-c(sum(datar$recruits[index])/sum(datar$stock[index]),mean(datar$recruits[c(length(datar$recruits)/2):length(datar$recruits)]))
     temp<-nls(datar$recruits~(A*datar$stock)/(1+(A*datar$stock)/B),start=list(A=parms[1],B=parms[2]))
     parms<-as.numeric(coef(temp))
    tres<-residuals(temp)
     ar1<-arima(tres,order=c(1,0,0),include.mean=FALSE)
     parms<-c(parms,as.numeric(ar1$coef))
   }
   if(select==2){
     if(is.null(initial$BHA)) stop("starting value for parameter BHA is missing")
     if(is.null(initial$BHB)) stop("starting value for parameter BHB is missing")
     if(is.null(initial$BHrho)) stop("starting value for parameter BHrho is missing")
      parms<-c(initial$BHA,initial$BHB,initial$BHrho)
   }
    bhnormcorr<-function(parms){
		 A<-parms[1]
     B<-parms[2]
     phi<-parms[3]    
     phi<-ifelse(phi>1,1,ifelse(phi<=-1,-1,phi))
     vpred<-(A*datar$stock)/(1+(A*datar$stock)/B) 
     res<-datar$recruits-vpred
     es<-c(res[1:c(length(res)-1)]*phi)
     sigma2w<-sum((res[-1]-es)^2)/length(res[-1])
     sump1<-0
     for(k in 2:n) sump1<-sump1+(datar$recruits[k]-phi*datar$recruits[k-1]-vpred[k]+phi*vpred[k-1])^2
     NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2w))-0.5*log(1-phi^2)+
             1/(2*sigma2w)*sump1+((1-phi^2)/(2*sigma2w))*(datar$recruits[1]-vpred[1])^2           
     NLL  
    }  
    model7 <- try(optim(parms,bhnormcorr,control=control),silent=TRUE)
 if(class(model7)!="try-error"){
      K<-length(model7$par)+1
      model7$hessian<-hessian(bhnormcorr,model7$par)
      model7$pred<-vpred<-(model7$par[1]*datar$stock)/(1+(model7$par[1]*datar$stock)/model7$par[2])
      res<-datar$recruits-vpred
      es<-c(res[1:c(length(res)-1)]*model7$par[3])
      model7$sigma2w<-sum((res[-1]-es)^2)/length(res[-1])
      model7$AICC<-2*model7$value+2*K*(n/(n-K-1))
      model7$resids<-datar$recruits-vpred
   }
}#BH normal correlated

############################### Model 8 #################################
if(any(model==8)){ # Beverton-Holt correlated log-normal errors
if(select==1){
  index<-which(datar$recruits>0)[1:2]
     parms<-c(sum(datar$recruits[index])/sum(datar$stock[index]),mean(datar$recruits[c(length(datar$recruits)/2):length(datar$recruits)]))
     tres<-log(datar$recruits)-log((parms[1]*datar$stock)/(1+((parms[1]*datar$stock)/parms[2])))
     ar1<-arima(tres,order=c(1,0,0),include.mean=FALSE)
     parms<-c(parms,as.numeric(ar1$coef))
   }
   if(select==2){
    if(is.null(initial$BHA)) stop("starting value for parameter BHA is missing")
    if(is.null(initial$BHB)) stop("starting value for parameter BHB is missing")
    if(is.null(initial$BHrho)) stop("starting value for parameter BHrho is missing")
    parms<-c(initial$BHA,initial$BHB,initial$BHrho)
   }
    bhlogcorr<-function(parms){
		  A<-parms[1]
      B<-parms[2]
      phi<-parms[3]
      phi<-ifelse(phi>1,1,ifelse(phi<=-1,-1,phi))
      vpred<-(A*datar$stock)/(1+(A*datar$stock)/B) 
      vpred<-ifelse(vpred<=0,0.0001,vpred)  
      res<-log(datar$recruits)-log(vpred)
      es<-c(res[1:c(length(res)-1)]*phi)
      sigma2w<-sum((res[-1]-es)^2)/length(res[-1])
      sump1<-0
      for(k in 2:n) sump1<-sump1+(log(datar$recruits[k])-phi*log(datar$recruits[k-1])-log(vpred[k])+phi*log(vpred[k-1])+(1-phi)*sigma2w/2)^2
         NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2w))+sum(log(datar$recruits))-0.5*log(1-phi^2)+
         1/(2*sigma2w)*sump1+((1-phi^2)/(2*sigma2w))*(log(datar$recruits[1])-log(vpred[1])+sigma2w/2)^2           
         NLL  
    }  
    model8 <- try(optim(parms,bhlogcorr,control=control),silent=TRUE)
   if(class(model8)!="try-error"){
      K<-length(model8$par)+1
      model8$hessian<-hessian(bhlogcorr,model8$par)
      model8$pred<-vpred<-(model8$par[1]*datar$stock)/(1+(model8$par[1]*datar$stock)/model8$par[2])
      res<-log(datar$recruits)-log(vpred)
      es<-c(res[1:c(length(res)-1)]*model8$par[3])
      model8$sigma2w<-sum((res[-1]-es)^2)/length(res[-1])  
      model8$AICC<-2*model8$value+2*K*(n/(n-K-1))
      model8$resids<-log(datar$recruits)-log(vpred)   
 }
  
}#BH log-normal correlated

########################## Model 9 ##################################
if(any(model==9)){ #Shepherd normal errors
   if(select==1){
     index<-which(datar$recruits>0)[1:2]
     parms<-c(sum(datar$recruits[index])/sum(datar$stock[index]),mean(datar$recruits)/2,1)
   }
   if(select==2){
     if(is.null(initial$SHA)) stop("starting value for parameter SHA is missing")
     if(is.null(initial$SHB)) stop("starting value for parameter SHB is missing")
     if(is.null(initial$SHC)) stop("starting value for parameter SHC is missing")
       parms<-c(initial$SHA,initial$SHB,initial$SHC)
   }
    shnorm<-function(parms){
		  A<-parms[1]
      B<-parms[2]
      C<-parms[3]
      vpred<-(A*datar$stock)/(1+((datar$stock/B)^C))        
      sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
      NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+
           1/(2*sigma2)*sum((datar$recruits-vpred)^2)
      NLL
      }  
    model9 <- try(optim(parms,shnorm,control=control),silent=FALSE)
if(class(model9)!="try-error"){
    K<-length(model9$par)+1
    model9$hessian<-hessian(shnorm,model9$par)
    model9$pred<-vpred<-(model9$par[1]*datar$stock)/(1+((datar$stock/model9$par[2])^model9$par[3]))
    model9$sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    model9$AICC<-2*model9$value+2*K*(n/(n-K-1))
    model9$resids<-datar$recruits-vpred
   }
}

################ Model 10 ##############################################
if(any(model==10)){ #Shepherd log-normal errors
 if(select==1){
   index<-which(datar$recruits>0)[1:2]
   parms<-c(sum(datar$recruits[index])/sum(datar$stock[index]),mean(datar$recruits)/2,1) 
   }
   if(select==2){
     if(is.null(initial$SHA)) stop("starting value for parameter SHA is missing")
     if(is.null(initial$SHB)) stop("starting value for parameter SHB is missing")
     if(is.null(initial$SHC)) stop("starting value for parameter SHC is missing")
     parms<-c(initial$SHA,initial$SHB,initial$SHC)
   }
    shlog<-function(parms){
      A<-parms[1]
      B<-parms[2]
      C<-parms[3]      
      vpred<-(A*datar$stock)/(1+((datar$stock/B)^C))
      vpred<-ifelse(vpred<=0,0.0001,vpred)
      sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
       NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+sum(log(datar$recruits))+
         1/(2*sigma2)*sum((log(datar$recruits)-log(vpred)+sigma2/2)^2)
      NLL      
         }  
    model10 <- try(optim(parms,shlog,control=control),silent=FALSE)
   if(class(model10)!="try-error"){
      K<-length(model10$par)+1
      model10$hessian<-hessian(shlog,model10$par)
      model10$pred<-vpred<-(model10$par[1]*datar$stock)/(1+((datar$stock/model10$par[2])^model10$par[3]))
      model10$sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
      model10$AICC<-2*model10$value+2*K*(n/(n-K-1))
      model10$resids<-log(datar$recruits)-log(vpred)
    }
  } # Shepherd log-normal
####################################Model 11 ####################################

if(any(model==11)){ #Deriso-SChnute normal errors
  if(select==1){
    index<-which(datar$recruits>0)[1:2]
    al<-sum(datar$recruits[index])/sum(datar$stock[index])
    newparms<-try(nls(datar$recruits~A*datar$stock*(1-B*C*datar$stock)^(1/C),
                      start=list(A=al,B=1/mean(datar$stock[c(length(datar$stock)/2):length(datar$stock)]),C=-0.5),
                      control=list(maxiter=10000)),silent=TRUE)
    if(class(newparms)!="try-error") parms<-c(as.numeric(summary(newparms)$coef[,1]))
    if(class(newparms)=="try-error") stop("Can't generate good starting values using automatic selection for Model 11. Try manual starting values.")
      }
  if(select==2){
    if(is.null(initial$DSA)) stop("starting value for parameter DSA is missing")
    if(is.null(initial$DSB)) stop("starting value for parameter DSB is missing")
    if(is.null(initial$DSC)) stop("starting value for parameter DSC is missing")
    parms<-c(initial$DSA,initial$DSB,initial$DSC)
  }
  dsnorm<-function(parms){
    A<-max(0,parms[1])
    B<-max(0,parms[2])
    C<-parms[3]
    vpred<-A*datar$stock*(1-B*C*datar$stock)^(1/C)  
    sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+
      1/(2*sigma2)*sum((datar$recruits-vpred)^2)
    NLL
  }  
   model11 <- try(optim(parms,dsnorm,control=control),silent=FALSE)
  if(class(model11)!="try-error"){
    K<-length(model11$par)+1
    model11$pred<-vpred<-(model11$par[1]*datar$stock)*(1-model11$par[2]*model11$par[3]*datar$stock)^(1/model11$par[3])        
    sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    model11$hessian<-hessian(dsnorm,model11$par)
    model11$sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    model11$AICC<-2*model11$value+2*K*(n/(n-K-1))
    model11$resids<-datar$recruits-vpred
  }
}

################ Model 12 ##############################################
if(any(model==12)){ #Deriso-Schnute log-normal errors
  if(select==1){
    index<-which(datar$recruits>0)[1:2]
    al<-sum(datar$recruits[index])/sum(datar$stock[index])
    newparms<-try(nls(datar$recruits~A*datar$stock*(1-B*C*datar$stock)^(1/C),
                      start=list(A=al,B=1/mean(datar$stock[c(length(datar$stock)/2):length(datar$stock)]),C=-0.5),
                      control=list(maxiter=10000)),silent=TRUE)
    if(class(newparms)!="try-error") parms<-c(as.numeric(summary(newparms)$coef[,1]))
    if(class(newparms)=="try-error") stop("Can't generate good starting values using automatic selection for Model 12. Try manual starting values.")
  }
  if(select==2){
    if(is.null(initial$DSA)) stop("starting value for parameter DSA is missing")
    if(is.null(initial$DSB)) stop("starting value for parameter DSB is missing")
    if(is.null(initial$DSC)) stop("starting value for parameter DSC is missing")
    parms<-c(initial$DSA,initial$DSB,initial$DSC)
  }
  dslnorm<-function(parms){
    A<-max(0,parms[1])
    B<-max(0,parms[2])
    C<-parms[3]
    vpred<-A*datar$stock*(1-B*C*datar$stock)^(1/C) 
    vpred<-ifelse(vpred<=0,0.0001,vpred)
    sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
    NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+sum(log(datar$recruits))+
      1/(2*sigma2)*sum((log(datar$recruits)-log(vpred)+sigma2/2)^2)
    NLL     
  }  
    model12 <- try(optim(parms,dslnorm,control=control),silent=FALSE)
  if(class(model12)!="try-error"){
    K<-length(model12$par)+1
    model12$hessian<-hessian(dslnorm,model12$par)
    model12$pred<-vpred<-model12$par[1]*datar$stock*(1-model12$par[2]*model12$par[3]*datar$stock)^(1/model12$par[3])
    model12$sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
    model12$AICC<-2*model12$value+2*K*(n/(n-K-1))
    model12$resids<-log(datar$recruits)-log(vpred)
  }
} # Deriso-Schnute log-normal

####################################Model 13 ####################################

if(any(model==13)){ #Myers depensation normal errors
  if(select==1){
    indexp<-which(datar$recruits>0)[1:2]
      al<-sum(datar$recruits[indexp])/sum(datar$stock[indexp])
      newparms<-try(nls(datar$recruits~(A*datar$stock^C)/(1+(datar$stock^C)/B),
                  start=list(A=al,B=mean(datar$recruits[c(length(datar$recruits)/2):length(datar$recruits)]),C=1),
                  control=list(maxiter=10000)),silent=TRUE)
    if(class(newparms)!="try-error") parms<-c(as.numeric(summary(newparms)$coef[1,1]),as.numeric(summary(newparms)$coef[2:3,1]))
    if(class(newparms)=="try-error") stop("Can't generate good starting values using automatic selection for Model 13. Try manual starting values.")
      
  }
  if(select==2){
    if(is.null(initial$MYA)) stop("starting value for parameter MYA is missing")
    if(is.null(initial$MYB)) stop("starting value for parameter MYB is missing")
    if(is.null(initial$MYC)) stop("starting value for parameter MYC is missing")
    parms<-c(initial$MYA,initial$MYB,initial$MYC)
  }
  mynorm<-function(parms){
    A<-max(0,parms[1])
    B<-max(0,parms[2])
    C<-max(0,parms[3])
    vpred<-(A*datar$stock^C)/(1+((datar$stock^C)/B))  
    sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+
      1/(2*sigma2)*sum((datar$recruits-vpred)^2)
    NLL
  }  
  model13 <- try(optim(parms,mynorm,control=control),silent=FALSE)
  if(class(model13)!="try-error"){
    K<-length(model13$par)+1
    model13$pred<-vpred<-(model13$par[1]*datar$stock^model13$par[3])/(1+((datar$stock^model13$par[3])/model13$par[2]))  
    sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    model13$hessian<-hessian(mynorm,model13$par)
    model13$sigma2<-sum((datar$recruits-vpred)^2)/length(datar$recruits)
    model13$AICC<-2*model13$value+2*K*(n/(n-K-1))
    model13$resids<-datar$recruits-vpred
  }
}

################ Model 14 ##############################################
if(any(model==14)){ #Myers Depensatory log-normal errors
  if(select==1){
    indexp<-which(datar$recruits>0)[1:2]
    al<-sum(datar$recruits[indexp])/sum(datar$stock[indexp])
    newparms<-try(nls(datar$recruits~(A*datar$stock^C)/(1+(datar$stock^C)/B),
                      start=list(A=al,B=mean(datar$recruits[c(length(datar$recruits)/2):length(datar$recruits)]),C=1),
                      control=list(maxiter=10000)),silent=TRUE)
    if(class(newparms)!="try-error") parms<-c(as.numeric(summary(newparms)$coef[1,1]),as.numeric(summary(newparms)$coef[2:3,1]))
    if(class(newparms)=="try-error") stop("Can't generate good starting values using automatic selection for Model 14. Try manual starting values.")
  }
  if(select==2){
    if(is.null(initial$MYA)) stop("starting value for parameter MYA is missing")
    if(is.null(initial$MYB)) stop("starting value for parameter MYB is missing")
    if(is.null(initial$MYC)) stop("starting value for parameter MYC is missing")
    parms<-c(initial$MYA,initial$MYB,initial$MYC)
  }
  mylnorm<-function(parms){
    A<-max(0,parms[1])
    B<-max(0,parms[2])
    C<-max(0,parms[3])
    vpred<-(A*datar$stock^C)/(1+((datar$stock^C)/B))  
    sigma2<-sum((log(datar$recruits)-log(vpred))^2)/length(datar$recruits)
    NLL<-n/2*log(2*pi)+n*log(sqrt(sigma2))+sum(log(datar$recruits))+
      1/(2*sigma2)*sum((log(datar$recruits)-log(vpred)+sigma2/2)^2)
    NLL     
  }  
  model14 <- try(optim(parms,mylnorm,control=control),silent=FALSE)
  if(class(model14)!="try-error"){
    K<-length(model14$par)+1
    model14$hessian<-hessian(mylnorm,model14$par)
    model14$pred<-vpred<-(model14$par[1]*datar$stock^model14$par[3])/(1+((datar$stock^model14$par[3])/model14$par[2]))  
    model14$sigma2<-sum((log(datar$recruits+0.0001)-log(vpred+0.0001))^2)/length(datar$recruits)
    model14$AICC<-2*model14$value+2*K*(n/(n-K-1))
    model14$resids<-log(datar$recruits+0.0001)-log(vpred+0.0001)
  }
} # Myers Depensation log-normal

######################### Consolidate Results
modelpt<-as.data.frame(matrix(0,ncol=c(length(model)),nrow=12))
rownames(modelpt)<-c("A","SE A","B","SE B","C","SE C","rho","SE rho","sigma2","sigma2w","NLL","AICC")
model_residuals<-as.data.frame(matrix(0,ncol=c(length(model)),nrow=length(datar$recruits)))
parm_correl<-as.data.frame(matrix(0,ncol=c(length(model)),
nrow=5))
rownames(parm_correl)<-c("A vs B","A vs C","A vs rho",
                         "B vs C", "B vs rho")
converge<-as.data.frame(matrix(-1,ncol=c(length(model)),nrow=1))
predictions<-as.data.frame(matrix(0,ncol=c(length(model)),nrow=length(datar$recruits)))

cnt<-0
modlabels<-rep("NA",length(model))

if(any(model==0)){
cnt<-cnt+1
modlabels[cnt]<-"Dens-Ind"
if(class(model0)!="try-error"){
 modelpt[c(1,2,9,11,12),cnt]<-c(model0$estimate,sqrt(solve(model0$hessian)),model0$sigma2,
   round(model0$minimum,2),model0$AICC)
   model_residuals[,cnt]<-model0$resids
   converge[cnt]<-model0$convergence
   predictions[,cnt]<-model0$pred
  }
  if(class(model0)=="try-error") modelpt[c(1,2,9,11,12),cnt]<-rep(NaN,7)
}
if(any(model==1)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"Ricker N-U"
  if(class(model1)!="try-error"){
    ses<-try(round(sqrt(diag(solve(model1$hessian))),10),silent=TRUE)
   if(class(ses)=="try-error") ses<-rep(NaN,length(model1$par))
     modelpt[c(1,2,3,4,9,11,12),cnt]<-c(model1$par[1],ses[1],model1$par[2],
      ses[2],model1$sigma2,round(model1$value,2),model1$AICC)
     model_residuals[,cnt]<-model1$resids
     corr<-try(cov2cor(solve(model1$hessian)),silent=TRUE)
     if(class(corr)!="try-error"){
       corr<-corr[lower.tri(corr)]
       parm_correl[1,cnt]<-corr[1]
     }
     
     if(class(corr)=="try-error"){
       parm_correl[1,cnt]<-NA
     }
     converge[cnt]<-model1$convergence
     predictions[,cnt]<-model1$pred 
  }
  if(class(model1)=="try-error") modelpt[c(1,2,3,4,7,11,12),cnt]<-rep(NaN,7)
}
if(any(model==2)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"Ricker L-U"
  if(class(model2)!="try-error"){
   ses<-try(sqrt(diag(solve(model2$hessian))),silent=TRUE)
   if(class(ses)=="try-error") ses<-rep(NaN,length(model2$par))
     modelpt[c(1,2,3,4,9,11,12),cnt]<-c(model2$par[1],ses[1],model2$par[2],
    ses[2],model2$sigma2,round(model2$value,2),model2$AICC)
    model_residuals[,cnt]<-model2$resids
     corr<-try(cov2cor(solve(model2$hessian)),silent=TRUE)
      if(class(corr)!="try-error"){
        corr<-corr[lower.tri(corr)]
        parm_correl[1,cnt]<-corr[1]
      }
      
      if(class(corr)=="try-error"){
        parm_correl[1,cnt]<-NA
      }
     converge[cnt]<-model2$convergence
     predictions[,cnt]<-model2$pred
   }
     if(class(model2)=="try-error") modelpt[c(1,2,3,4,7,11,12),cnt]<-rep(NaN,7)
}
if(any(model==3)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"Ricker N-C"
 if(class(model3)!="try-error"){
   ses<-try(sqrt(diag(solve(model3$hessian))),silent=TRUE)
   if(class(ses)=="try-error") ses<-rep(NaN,length(model3$par))
    modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-c(model3$par[1],ses[1],model3$par[2],
    ses[2],model3$par[3],ses[3],model3$sigma2w,round(model3$value,2),model3$AICC)
    model_residuals[,cnt]<-model3$resids
     corr<-try(cov2cor(solve(model3$hessian)),silent=TRUE)
     if(class(corr)!="try-error"){
       corr<-corr[lower.tri(corr)]
       parm_correl[1,cnt]<-corr[1]
       parm_correl[3,cnt]<-corr[2]
       parm_correl[5,cnt]<-corr[3]
     }
     if(class(corr)=="try-error"){
       parm_correl[1,cnt]<-NA
       parm_correl[3,cnt]<-NA
       parm_correl[5,cnt]<-NA
     }
    converge[cnt]<-model3$convergence
    predictions[,cnt]<-model3$pred
   }
  if(class(model3)=="try-error") modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-rep(NaN,9)
}
if(any(model==4)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"Ricker L-C"
 if(class(model4)!="try-error"){
   ses<-try(sqrt(diag(solve(model4$hessian))),silent=TRUE)
    if(class(ses)=="try-error") ses<-rep(NaN,length(model4$par))
   modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-c(model4$par[1],ses[1],model4$par[2],
    ses[2],model4$par[3],ses[3],model4$sigma2w,round(model4$value,2),model4$AICC)
   model_residuals[,cnt]<-model4$resids
   corr<-try(cov2cor(solve(model4$hessian)),silent=TRUE)
   if(class(corr)!="try-error"){
     corr<-corr[lower.tri(corr)]
     parm_correl[1,cnt]<-corr[1]
     parm_correl[3,cnt]<-corr[2]
     parm_correl[5,cnt]<-corr[3]
   }
   if(class(corr)=="try-error"){
     parm_correl[1,cnt]<-NA
     parm_correl[3,cnt]<-NA
     parm_correl[5,cnt]<-NA
    
   }
   converge[cnt]<-model4$convergence
   predictions[,cnt]<-model4$pred
  }
  if(class(model4)=="try-error") modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-rep(NaN,9)
}
if(any(model==5)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"BH N-U"
  if(class(model5)!="try-error"){
    ses<-try(sqrt(diag(solve(model5$hessian))),silent=TRUE)
   if(class(ses)=="try-error") ses<-rep(NaN,length(model5$par))
     modelpt[c(1,2,3,4,9,11,12),cnt]<-c(model5$par[1],ses[1],model5$par[2],
      ses[2],model5$sigma2,round(model5$value,2),model5$AICC)
    model_residuals[,cnt]<-model5$resids
 corr<-try(cov2cor(solve(model5$hessian)),silent=TRUE)
 if(class(corr)!="try-error"){
   corr<-corr[lower.tri(corr)]
   parm_correl[1,cnt]<-corr[1]
 }
 
 if(class(corr)=="try-error"){
   parm_correl[1,cnt]<-NA
 }
     converge[cnt]<-model5$convergence
     predictions[,cnt]<-model5$pred
  }
  if(class(model5)=="try-error") modelpt[c(1,2,3,4,7,11,12),cnt]<-rep(NaN,7)
}
if(any(model==6)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"BH L-U"
  if(class(model6)!="try-error"){
   ses<-try(sqrt(diag(solve(model6$hessian))),silent=TRUE)
   if(class(ses)=="try-error") ses<-rep(NaN,length(model6$par))
   modelpt[c(1,2,3,4,9,11,12),cnt]<-c(model6$par[1],ses[1],model6$par[2],
    ses[2],model6$sigma2,round(model6$value,2),model6$AICC)
   model_residuals[,cnt]<-model6$resids
     corr<-try(cov2cor(solve(model6$hessian)),silent=TRUE)
     if(class(corr)!="try-error"){
       corr<-corr[lower.tri(corr)]
       parm_correl[1,cnt]<-corr[1]
     }
     if(class(corr)=="try-error"){
       parm_correl[1,cnt]<-NA
     }
  converge[cnt]<-model6$convergence
  predictions[,cnt]<-model6$pred
  }
  if(class(model6)=="try-error") modelpt[c(1,2,3,4,7,11,12),cnt]<-rep(NaN,7)
}
if(any(model==7)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"BH N-C"
 if(class(model7)!="try-error"){
   ses<-try(sqrt(diag(solve(model7$hessian))),silent=TRUE)
   if(class(ses)=="try-error") ses<-rep(NaN,length(model7$par))
   modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-c(model7$par[1],ses[1],model7$par[2],
    ses[2],model7$par[3],ses[3],model7$sigma2w,round(model7$value,2),model7$AICC)
    model_residuals[,cnt]<-model7$resids
 corr<-try(cov2cor(solve(model7$hessian)),silent=TRUE)
 if(class(corr)!="try-error"){
   corr<-corr[lower.tri(corr)]
   parm_correl[1,cnt]<-corr[1]
   parm_correl[3,cnt]<-corr[2]
   parm_correl[5,cnt]<-corr[3]
 }
 if(class(corr)=="try-error"){
   parm_correl[1,cnt]<-NA
   parm_correl[3,cnt]<-NA
   parm_correl[5,cnt]<-NA
 }
   converge[cnt]<-model7$convergence
   predictions[,cnt]<-model7$pred
  }
  if(class(model7)=="try-error") modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-rep(NaN,9)
}
if(any(model==8)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"BH L-C"
 if(class(model8)!="try-error"){
   ses<-try(sqrt(diag(solve(model8$hessian))),silent=TRUE)
 if(class(ses)=="try-error") ses<-rep(NaN,length(model8$par))

   modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-c(model8$par[1],ses[1],model8$par[2],
    ses[2],model8$par[3],ses[3],model8$sigma2w,round(model8$value,2),model8$AICC)
   model_residuals[,cnt]<-model8$resids
 corr<-try(cov2cor(solve(model8$hessian)),silent=TRUE)
 if(class(corr)!="try-error"){
   corr<-corr[lower.tri(corr)]
   parm_correl[1,cnt]<-corr[1]
   parm_correl[3,cnt]<-corr[2]
   parm_correl[5,cnt]<-corr[3]
 }
 if(class(corr)=="try-error"){
   parm_correl[1,cnt]<-NA
   parm_correl[3,cnt]<-NA
   parm_correl[5,cnt]<-NA
 
 }
 converge[cnt]<-model8$convergence
predictions[,cnt]<-model8$pred
  }
  if(class(model8)=="try-error") modelpt[c(1,2,3,4,7,8,10,11,12),cnt]<-rep(NaN,9)
}
if(any(model==9)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"SH N-U"
 if(class(model9)!="try-error"){
   ses<-try(sqrt(diag(solve(model9$hessian))),silent=TRUE)
 if(class(ses)=="try-error") ses<-rep(NaN,length(model9$par))
   modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-c(model9$par[1],ses[1],model9$par[2],
    ses[2],model9$par[3],ses[3],model9$sigma2,round(model9$value,2),model9$AICC)
   model_residuals[,cnt]<-model9$resids
   corr<-try(cov2cor(solve(model9$hessian)),silent=TRUE)
   if(class(corr)!="try-error"){
     corr<-corr[lower.tri(corr)]
     parm_correl[1,cnt]<-corr[1]
     parm_correl[2,cnt]<-corr[2]
     parm_correl[4,cnt]<-corr[3]
   }
   if(class(corr)=="try-error"){
     parm_correl[1,cnt]<-NA
     parm_correl[2,cnt]<-NA
     parm_correl[4,cnt]<-NA
   }
   converge[cnt]<-model9$convergence
   predictions[,cnt]<-model9$pred
  }
  if(class(model9)=="try-error") modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-rep(NaN,9)
}

if(any(model==10)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"SH L-U"
 if(class(model10)!="try-error"){
   ses<-try(sqrt(diag(solve(model10$hessian))),silent=TRUE)
 if(class(ses)=="try-error") ses<-rep(NaN,length(model10$par))
   modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-c(model10$par[1],ses[1],model10$par[2],
    ses[2],model10$par[3],ses[3],model10$sigma2,round(model10$value,2),model10$AICC)
   model_residuals[,cnt]<-model10$resids
   corr<-try(cov2cor(solve(model10$hessian)),silent=TRUE)
   if(class(corr)!="try-error"){
     corr<-corr[lower.tri(corr)]
     parm_correl[1,cnt]<-corr[1]
     parm_correl[2,cnt]<-corr[2]
     parm_correl[4,cnt]<-corr[3]
   
   }
   if(class(corr)=="try-error"){
     parm_correl[1,cnt]<-NA
     parm_correl[2,cnt]<-NA
     parm_correl[4,cnt]<-NA
   }
    converge[cnt]<-model10$convergence
    predictions[,cnt]<-model10$pred
  }
  if(class(model10)=="try-error") modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-rep(NaN,9)
}

if(any(model==11)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"DS N-U"
  if(class(model11)!="try-error"){
    ses<-try(sqrt(diag(solve(model11$hessian))),silent=TRUE)
    if(class(ses)=="try-error") ses<-rep(NaN,length(model11$par))
    modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-c(model11$par[1],ses[1],model11$par[2],
                                         ses[2],model11$par[3],ses[3],model11$sigma2,round(model11$value,2),model11$AICC)
    model_residuals[,cnt]<-model11$resids
    corr<-try(cov2cor(solve(model11$hessian)),silent=TRUE)
    if(class(corr)!="try-error"){
      corr<-corr[lower.tri(corr)]
      parm_correl[1,cnt]<-corr[1]
      parm_correl[2,cnt]<-corr[2]
      parm_correl[4,cnt]<-corr[3]
      
    }
    if(class(corr)=="try-error"){
      parm_correl[1,cnt]<-NA
      parm_correl[2,cnt]<-NA
      parm_correl[4,cnt]<-NA
    }
    converge[cnt]<-model11$convergence
    predictions[,cnt]<-model11$pred
  }
  if(class(model11)=="try-error") modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-rep(NaN,9)
}


if(any(model==12)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"DS L-U"
  if(class(model12)!="try-error"){
    ses<-try(sqrt(diag(solve(model12$hessian))),silent=TRUE)
    if(class(ses)=="try-error") ses<-rep(NaN,length(model12$par))
    modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-c(model12$par[1],ses[1],model12$par[2],
                                         ses[2],model12$par[3],ses[3],model12$sigma2,round(model12$value,2),model12$AICC)
    model_residuals[,cnt]<-model12$resids
    corr<-try(cov2cor(solve(model12$hessian)),silent=TRUE)
    if(class(corr)!="try-error"){
      corr<-corr[lower.tri(corr)]
      parm_correl[1,cnt]<-corr[1]
      parm_correl[2,cnt]<-corr[2]
      parm_correl[4,cnt]<-corr[3]
      
    }
    if(class(corr)=="try-error"){
      parm_correl[1,cnt]<-NA
      parm_correl[2,cnt]<-NA
      parm_correl[4,cnt]<-NA
    }
    converge[cnt]<-model12$convergence
    predictions[,cnt]<-model12$pred
  }
  if(class(model12)=="try-error") modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-rep(NaN,9)
}


if(any(model==13)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"MY N-U"
  if(class(model13)!="try-error"){
    ses<-try(sqrt(diag(solve(model13$hessian))),silent=TRUE)
    if(class(ses)=="try-error") ses<-rep(NaN,length(model13$par))
    modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-c(model13$par[1],ses[1],model13$par[2],
                                         ses[2],model13$par[3],ses[3],model13$sigma2,round(model13$value,2),model13$AICC)
    model_residuals[,cnt]<-model13$resids
    corr<-try(cov2cor(solve(model13$hessian)),silent=TRUE)
    if(class(corr)!="try-error"){
      corr<-corr[lower.tri(corr)]
      parm_correl[1,cnt]<-corr[1]
      parm_correl[2,cnt]<-corr[2]
      parm_correl[4,cnt]<-corr[3]
      
    }
    if(class(corr)=="try-error"){
      parm_correl[1,cnt]<-NA
      parm_correl[2,cnt]<-NA
      parm_correl[4,cnt]<-NA
    }
    converge[cnt]<-model13$convergence
    predictions[,cnt]<-model13$pred
  }
  if(class(model13)=="try-error") modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-rep(NaN,9)
}

if(any(model==14)){ 
  cnt<-cnt+1
  modlabels[cnt]<-"MY L-U"
  if(class(model14)!="try-error"){
    ses<-try(sqrt(diag(solve(model14$hessian))),silent=TRUE)
    if(class(ses)=="try-error") ses<-rep(NaN,length(model14$par))
    modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-c(model14$par[1],ses[1],model14$par[2],
                                         ses[2],model14$par[3],ses[3],model14$sigma2,round(model14$value,2),model14$AICC)
    model_residuals[,cnt]<-model14$resids
    corr<-try(cov2cor(solve(model14$hessian)),silent=TRUE)
    if(class(corr)!="try-error"){
      corr<-corr[lower.tri(corr)]
      parm_correl[1,cnt]<-corr[1]
      parm_correl[2,cnt]<-corr[2]
      parm_correl[4,cnt]<-corr[3]
      
    }
    if(class(corr)=="try-error"){
      parm_correl[1,cnt]<-NA
      parm_correl[2,cnt]<-NA
      parm_correl[4,cnt]<-NA
    }
    converge[cnt]<-model14$convergence
    predictions[,cnt]<-model14$pred
  }
  if(class(model14)=="try-error") modelpt[c(1,2,3,4,5,6,9,11,12),cnt]<-rep(NaN,9)
}

colnames(modelpt)<-colnames(model_residuals)<-colnames(parm_correl)<-colnames(converge)<-colnames(predictions)<-modlabels
evidence_ratios<-as.data.frame(matrix(0,nrow=c(length(model)),ncol=4))
colnames(evidence_ratios)<-c("Model","Di","Akaike Wgt","Evidence_Ratio")
evidence_ratios[,1]<-modlabels
evidence_ratios[,2]<-as.numeric(modelpt[12,]-min(modelpt[12,]))
delt<-as.numeric(exp(-0.5*(evidence_ratios[,2])))
evidence_ratios[,3]<-delt/sum(delt)
ratios<-ifelse(evidence_ratios[,3]<=0,0,max(evidence_ratios[,3])/evidence_ratios[,3])
evidence_ratios[,4]<-ratios
evidence_ratios<-evidence_ratios[order(evidence_ratios[,3],decreasing=TRUE),]
evidence_ratios[,2]<-round(evidence_ratios[,2],3)
evidence_ratios[,3]<-round(evidence_ratios[,3],3)
evidence_ratios[,4]<-round(evidence_ratios[,4],3)
if(plot==TRUE){
plot(datar$recruits~datar$stock,ylim=c(0,max(predictions)*1.20),
xlim=c(0,max(datar$stock)*1.20),type="p",
xlab="Spawning Stock (or Eggs)",ylab="Recruits")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
"gray80")
cols<-rainbow(length(model))
ltys<-c(sample(1:6,length(model),replace=TRUE))
for(cc in 1:ncol(predictions)){
  lines(predictions[,cc]~datar$stock,col=cols[cc],lty=ltys[cc],lwd=2)
}
points(datar$recruits~datar$stock,pch=16)
legend("bottomright",legend=colnames(predictions),bty="n",
col=cols,cex=0.6,lty=ltys,
lwd=2)
}
results<-list(results=modelpt,evidence_ratios=evidence_ratios,convergence=converge,correlations=parm_correl,
              predicted=predictions,residuals=model_residuals)
return(results)
}#function end



