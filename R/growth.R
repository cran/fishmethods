growth<-function(intype=1,unit=1,size=NULL,age=NULL,calctype=1,wgtby=1,se2=NULL,error=1,
                 specwgt=0.0001,
         Sinf=NULL,K=NULL,t0=NULL,B=3,graph=TRUE,
         control=list(maxiter=10000,minFactor=1/1024,tol=1e-5)){
   if(is.null(size)) 
         stop ("size is missing") 
   if(is.null(age)) 
         stop ("age is missing") 
   if(is.null(unit)) 
         stop ("unit is missing.") 
   if(length(age)!=length(size)) stop ("Vectors of different lengths")
   if((is.null(Sinf)|is.null(K)|is.null(t0))) stop("Values of Sinf, K, and t0 must be specified")
   if(intype==2 & wgtby==2 & is.null(se2)) stop("se2: You need to specify the se2 vector")
   if(intype==1 & calctype==1 & wgtby==2) stop("wgtby: You can't weight an individual observation, only means")
   if(intype==2 & calctype==2) stop("calctype: Data already inputted as means")
   wgts<-NULL;d4<-NULL;x<-NULL;vbl<-NULL;gomp<-NULL;logist<-NULL;vpred<-NULL;tempred<-NULL
   if(intype==1){
		if(calctype==1){
			  x<-as.data.frame(cbind(size,age)) 
   			  x<-x[!is.na(x$size) & !is.na(x$age),]
                wgts<-rep(1,length(x$size))
            }
		if(calctype==2){
              x<-as.data.frame(cbind(size,age)) 
		    x<-x[!is.na(x$size) & !is.na(x$age),]
		    d4<-merge(aggregate(x$size,list(x$age),mean),
                aggregate(x$size,list(x$age),function(x){var(x)/length(x)}),by.y=c("Group.1"),
                by.x=c("Group.1"))
                names(d4)<-c("age","size","var")
              x<-d4[!is.na(d4$size) & !is.na(d4$age),]
              x$var<-ifelse(x$var<=0.0|is.na(x$var),NA,x$var)
              if(wgtby==2) wgts<-1/x$var
              if(wgtby==1) wgts<-rep(1,length(x$var))
              wgts<-ifelse(is.na(wgts),specwgt,wgts)
            }
    }

   if(intype==2){
		if(wgtby==1){
	           x<-as.data.frame(cbind(size,age)) 
   		      x<-x[!is.na(x$size) & !is.na(x$age),]
                wgts<-rep(1,length(x$size))
           }
         if(wgtby==2){
             x<-as.data.frame(cbind(size,age,se2)) 
   		        x<-x[!is.na(x$size) & !is.na(x$age)& !is.na(x$se2),]
             wgts<-1/x$se2
     	}
     }

if(unit==1){#Length
   	   if(error==1){#additive
  	 	vbl<-try(nls(size~Sinf*(1-exp(-(K*(age-t0)))),data=x,       
        	weights=wgts,start=list(Sinf=Sinf,K=K,t0=t0),control=control),silent=FALSE)
  	 	gomp<-try(nls(size~Sinf*exp(-exp(-K*(age-t0))),data=x,       
  	 	              weights=wgts,start=list(Sinf=Sinf,K=K,t0=t0),control=control),silent=FALSE)
  	 	logist<-try(nls(size~Sinf/(1+exp(-K*(age-t0))),data=x,       
  	 	                weights=wgts,start=list(Sinf=Sinf,K=K,t0=t0),control=control),silent=FALSE)
  	 	
        }
	   if(error==2){#multiplicative
 		      mult1<-function(size,age,Sinf,K,t0){
        			vpred<-Sinf*(1-exp(-K*(age-t0)))
              		vpred<-ifelse(vpred<0,abs(vpred),vpred)
        	 		sqrt(wgts)*(log(size)-log(vpred)) 
      		}
    		vbl <- try(nls( ~ mult1(size,age,Sinf,K,t0), data = x, 
                  start = list(Sinf = Sinf, K = K,t0=t0),control=control, 
                   trace = FALSE),silent=FALSE) 
    		
    		mult3<-function(size,age,Sinf,K,t0){
    		  vpred<-Sinf*exp(-exp(-K*(age-t0)))
    		  vpred<-ifelse(vpred<0,abs(vpred),vpred)
    		  sqrt(wgts)*(log(size)-log(vpred)) 
    		}
    		gomp <- try(nls( ~ mult3(size,age,Sinf,K,t0), data = x, 
    		                 start = list(Sinf = Sinf, K = K,t0=t0),,control=control, 
    		                 trace = FALSE),silent=FALSE) 
    		# Logistic
    		mult4<-function(size,age,Sinf,K,t0){
    		  vpred<-Sinf/(1+exp(-K*(age-t0)))
    		  vpred<-ifelse(vpred<0,abs(vpred),vpred)
    		  sqrt(wgts)*(log(size)-log(vpred)) 
    		}
    		logist <- try(nls( ~ mult4(size,age,Sinf,K,t0), data = x, 
    		                   start = list(Sinf = Sinf, K = K,t0=t0),control=control,
    		                   trace = FALSE),silent=FALSE) 
       }        
 }#unit==1
       

if(unit==2){#Weight
   	   if(error==1){#additive
           x$Bs<-rep(B,times=length(x$size))
  	 	vbl<-try(nls(size~Sinf*((1-exp(-K*(age-t0)))^Bs),data=x,       
        		weights=wgts,start=list(Sinf=Sinf,K=K,t0=t0),control=control),silent=FALSE)
             
		  gomp<-try(nls(size~Sinf*exp(-exp(-K*(age-t0))),data=x,       
        	     weights=wgts,start=list(Sinf=Sinf,K=K,t0=t0),control=control),silent=FALSE)
             
          logist<-try(nls(size~Sinf/(1+exp(-K*(age-t0))),data=x,       
        		weights=wgts,start=list(Sinf=Sinf,K=K,t0=t0),control=control),silent=FALSE)

        }#error=1
	   if(error==2){#multiplicative
		# Von Bertalanffy
           x$Bs<-rep(B,times=length(x$size))
		 mult2<-function(size,age,Bs,Sinf,K,t0){
        			vpred<-Sinf*((1-exp(-K*(age-t0)))^Bs)
              		vpred<-ifelse(vpred<0,abs(vpred),vpred)
        	 		sqrt(wgts)*(log(size)-log(vpred)) 
      		}
    		 vbl <- try(nls( ~ mult2(size,age,Bs,Sinf,K,t0), data = x, 
                start = list(Sinf = Sinf, K = K,t0=t0),control=control, 
                trace = FALSE),silent=FALSE)
 		 
		# Gompertz
 		mult3<-function(size,age,Sinf,K,t0){
        			vpred<-Sinf*exp(-exp(-K*(age-t0)))
              		vpred<-ifelse(vpred<0,abs(vpred),vpred)
        	 		sqrt(wgts)*(log(size)-log(vpred)) 
      		}
    		gomp <- try(nls( ~ mult3(size,age,Sinf,K,t0), data = x, 
                 start = list(Sinf = Sinf, K = K,t0=t0),control=control,
                 trace = FALSE),silent=FALSE) 
 		
		# Logistic
 		mult4<-function(size,age,Sinf,K,t0){
        			vpred<-Sinf/(1+exp(-K*(age-t0)))
              		vpred<-ifelse(vpred<0,abs(vpred),vpred)
        	 		sqrt(wgts)*(log(size)-log(vpred)) 
      		}
    		logist <- try(nls( ~ mult4(size,age,Sinf,K,t0), data = x, 
                   start = list(Sinf = Sinf, K = K,t0=t0),control=control,
                   trace = FALSE),silent=FALSE) 
          }#error=2
 }  #unit=2


if(unit==1){
  if(intype==1){
    if(calctype==1){
      if(unit==1) lab<-c("Response: Individual Length")
      if(unit==2) lab<-c("Response: Individual Weight")
      lab<-c(lab,"Least Squares: Unweighted")
    }
    if(calctype==2) {
      if(unit==1) lab<-c("Response: Mean Length")
      if(unit==2) lab<-c("Response: Mean Weight")
      if(wgtby==1) lab<-c(lab,"Least Squares: Unweighted")
      if(wgtby==2) lab<-c(lab,"Least Squares: Weighted")
    }
  }
  if(intype==2){
    if(unit==1) lab<-c("Response: Mean Length")
    if(unit==2) lab<-c("Response: Mean Weight")
    if(wgtby==1) lab<-c(lab,"Least Squares: Unweighted")
    if(wgtby==2) lab<-c(lab,"Least Squares: Weighted")
  }
  if(error==1){
    lab1<-c("Model: Sinf*(1-exp(-K*(t-t0)))+e")
    lab2<-c("Model: Sinf*exp(-exp(-K*(age-t0)))+e")
    lab3<-c("Model: Sinf/(1+exp(-K*(age-t0))+e")
  }
  if(error==2){
    lab1<-c("Model: Sinf*(1-exp(-K*(t-t0)))*exp(e)")
    lab2<-c("Model: Sinf*exp(-exp(-K*(age-t0)))*exp(e)")
    lab3<-c("Model: Sinf/(1+exp(-K*(age-t0)))*exp(e)")
  }
  
  if(graph==TRUE){
    par(mfrow=c(3,2))
    ylab1<-ifelse(unit==1,"Length","Weight")
    ylab2<-ifelse(error==1,"Residual","Residual (log)")
    
    #Plot graphs
    if(class(vbl)!="try-error"){ 
      plot(x$size~x$age,col="blue",ylab=ylab1,xlab="Age",main="Von Bertalanffy Fit",xlim=c(0,max(x$age)+1))
      lines(summary(vbl)$parameters[1]*((1-exp(-summary(vbl)$parameters[2]*(
        seq(0,max(x$age)+1,1)-summary(vbl)$parameters[3]))))~seq(0,max(x$age)+1,1),
        col="red",lwd=2)  
      if(error==1){
        plot(summary(vbl)$residuals~x$age,ylim=c(-max(abs(summary(vbl)$residuals)),
                                                 max(abs(summary(vbl)$residuals))),col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
             main="Von Bertalanffy Residuals")
      }
      if(error==2){
        tempred<-summary(vbl)$parameters[1]*((1-exp(-summary(vbl)$parameters[2]*(
          x$age-summary(vbl)$parameters[3]))))
        plot(c(log(x$size)-log(tempred))~x$age,col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
             main="Von Bertalanffy Residuals")
      }
      abline(h=0,col="red",lwd=2)
    }
    
    if(class(gomp)!="try-error"){ 
      plot(x$size~x$age,col="blue",ylab=ylab1,xlab="Age",main="Gompertz Fit",xlim=c(0,max(x$age)+1))
      lines(summary(gomp)$parameters[1]*exp(-exp(-summary(gomp)$parameters[2]*(
        seq(0,max(x$age)+1,1)-summary(gomp)$parameters[3])))~seq(0,max(x$age)+1,1),
        col="red",lwd=2)  
      if(error==1){
        plot(summary(gomp)$residuals~x$age,ylim=c(-max(abs(summary(gomp)$residuals)),
                                                  max(abs(summary(gomp)$residuals))),col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
             main="Gompertz Residuals")
      }
      if(error==2){
        tempred<-summary(gomp)$parameters[1]*exp(-exp(-summary(gomp)$parameters[2]*(
          x$age-summary(gomp)$parameters[3])))
        plot(c(log(x$size)-log(tempred))~x$age,col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
             main="Gompertz Residuals")
      }
      abline(h=0,col="red",lwd=2)
    }
    
    if(class(logist)!="try-error"){ 
      plot(x$size~x$age,col="blue",ylab=ylab1,xlab="Age",main="Logistic Fit",xlim=c(0,max(x$age)+1))
      lines(summary(logist)$parameters[1]/(1+exp(-summary(logist)$parameters[2]*(
        seq(0,max(x$age)+1,1)-summary(logist)$parameters[3])))~seq(0,max(x$age)+1,1),
        col="red",lwd=2)  
      if(error==1){
        plot(summary(logist)$residuals~x$age,ylim=c(-max(abs(summary(logist)$residuals)),
                                                    max(abs(summary(logist)$residuals))),col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
             main="Logisitic Residuals")
      }
      if(error==2){
        tempred<-summary(logist)$parameters[1]/(1+exp(-summary(logist)$parameters[2]*(
          x$age-summary(logist)$parameters[3])))
        plot(c(log(x$size)-log(tempred))~x$age,col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
             main="Logistic Residuals")
      }
      abline(h=0,col="red",lwd=2)
    }
  }#graph true  
  
  if(class(vbl)=="try-error") vbl<-"Fit failed"
  if(class(gomp)=="try-error") gomp<-"Fit failed"
  if(class(logist)=="try-error") logist<-"Fit failed"
  nlsout<-list(c(lab1,lab),vbl,c(lab2,lab),gomp,c(lab3,lab),logist)
  names(nlsout)<-c("vonbert","vout","gompertz","gout","logistic","lout")
  par(mfrow=c(1,1))
  return(nlsout)
 }
if(unit==2){
      if(intype==1){
      	if(calctype==1){
       	   if(unit==1) lab<-c("Response: Individual Length")
       	   if(unit==2) lab<-c("Response: Individual Weight")
		   lab<-c(lab,"Least Squares: Unweighted")
       	}
     	 if(calctype==2) {
		   if(unit==1) lab<-c("Response: Mean Length")
		   if(unit==2) lab<-c("Response: Mean Weight")
       	   if(wgtby==1) lab<-c(lab,"Least Squares: Unweighted")
       	   if(wgtby==2) lab<-c(lab,"Least Squares: Weighted")
    		 }
      }
	if(intype==2){
	 	   if(unit==1) lab<-c("Response: Mean Length")
		   if(unit==2) lab<-c("Response: Mean Weight")
       	   if(wgtby==1) lab<-c(lab,"Least Squares: Unweighted")
       	   if(wgtby==2) lab<-c(lab,"Least Squares: Weighted")
    	 }
    	 if(error==1){
      
          lab1<-c("Model: Sinf*(1-exp(-K*(t-t0)))^B+e")
		 lab2<-c("Model: Sinf*exp(-exp(-K*(age-t0))+e")
		 lab3<-c("Model: Sinf/(1+exp(-K*(age-t0))+e")
        }
	 if(error==2){
         lab1<-c("Model: (Sinf*(1-exp(-K*(t-t0)))^B)*exp(e)")
         lab2<-c("Model: Sinf*exp(-exp(-K*(age-t0)))*exp(e)")
         lab3<-c("Model: Sinf/(1+exp(-K*(age-t0)))*exp(e)")
       }

 if(graph==TRUE){
       par(mfrow=c(3,2))
       ylab1<-ifelse(unit==1,"Length","Weight")
       ylab2<-ifelse(error==1,"Residual","Residual (log)")

    		#Plot graphs
    		if(class(vbl)!="try-error"){ 
      		plot(x$size~x$age,col="blue",ylab=ylab1,xlab="Age",main="Von Bertalanffy Fit",xlim=c(0,max(x$age)+1))
	 		lines(summary(vbl)$parameters[1]*((1-exp(-summary(vbl)$parameters[2]*(
       			seq(0,max(x$age)+1,1)-summary(vbl)$parameters[3]))))^B~seq(0,max(x$age)+1,1),
      			 col="red",lwd=2)  
       		if(error==1){
      			plot(summary(vbl)$residuals~x$age,ylim=c(-max(abs(summary(vbl)$residuals)),
          		max(abs(summary(vbl)$residuals))),col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
          		main="Von Bertalanffy Residuals")
       		}
       		if(error==2){
          		tempred<-summary(vbl)$parameters[1]*((1-exp(-summary(vbl)$parameters[2]*(
          		x$age-summary(vbl)$parameters[3])))^B)
      			plot(c(log(x$size)-log(tempred))~x$age,col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
          		main="Von Bertalanffy Residuals")
       		}
	   		abline(h=0,col="red",lwd=2)
         }

        if(class(gomp)!="try-error"){ 
      		plot(x$size~x$age,col="blue",ylab=ylab1,xlab="Age",main="Gompertz Fit",xlim=c(0,max(x$age)+1))
	 		lines(summary(gomp)$parameters[1]*exp(-exp(-summary(gomp)$parameters[2]*(
       		seq(0,max(x$age)+1,1)-summary(gomp)$parameters[3])))~seq(0,max(x$age)+1,1),
       		col="red",lwd=2)  
     		if(error==1){
      			plot(summary(gomp)$residuals~x$age,ylim=c(-max(abs(summary(gomp)$residuals)),
          		max(abs(summary(gomp)$residuals))),col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
          		main="Gompertz Residuals")
       		}
       		if(error==2){
          		tempred<-summary(gomp)$parameters[1]*exp(-exp(-summary(gomp)$parameters[2]*(
       				x$age-summary(gomp)$parameters[3])))
      			plot(c(log(x$size)-log(tempred))~x$age,col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
          		main="Gompertz Residuals")
       		}
	 		abline(h=0,col="red",lwd=2)
        }
     
	   if(class(logist)!="try-error"){ 
      		plot(x$size~x$age,col="blue",ylab=ylab1,xlab="Age",main="Logistic Fit",xlim=c(0,max(x$age)+1))
	 		lines(summary(logist)$parameters[1]/(1+exp(-summary(logist)$parameters[2]*(
       			seq(0,max(x$age)+1,1)-summary(logist)$parameters[3])))~seq(0,max(x$age)+1,1),
       			col="red",lwd=2)  
       		if(error==1){
      			plot(summary(logist)$residuals~x$age,ylim=c(-max(abs(summary(logist)$residuals)),
          		max(abs(summary(logist)$residuals))),col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
          		main="Logisitic Residuals")
        		}
        		if(error==2){
          		tempred<-summary(logist)$parameters[1]/(1+exp(-summary(logist)$parameters[2]*(
       				x$age-summary(logist)$parameters[3])))
      			plot(c(log(x$size)-log(tempred))~x$age,col="blue",ylab=ylab2,xlab="Age",xlim=c(0,max(x$age)+1),
          		main="Logistic Residuals")
        		}
	   		abline(h=0,col="red",lwd=2)
       }

    }#graph true  
      if(class(vbl)=="try-error") vbl<-"Fit failed"
 	 if(class(gomp)=="try-error") gomp<-"Fit failed"
      if(class(logist)=="try-error") logist<-"Fit failed"
	 nlsout<-list(c(lab1,lab),vbl,c(lab2,lab),gomp,c(lab3,lab),logist)
      names(nlsout)<-c("vonbert","vout","gompertz","gout","logistic","lout")
par(mfrow=c(1,1))
      return(nlsout)
 }
}# end of function
