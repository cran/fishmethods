convmort<-function(value=NULL,fromto=1,type=2,M=NULL){
  if(is.null(value)) stop("value is missing.")
  if(!any(fromto %in% c(1,2))) stop("'fromto' value incorrect.")
  if(!any(type %in% c(1,2))) stop("'type' value incorrect.")
  if(type==2 & is.null(M)) stop("Need M for a Type 2 fishery.")

  if(fromto==1){
    if(type==1){
      mu<-1-exp(-value)
      return(mu)
    }
    if(type==2){
      mu<-(value*(1-exp(-value-M)))/(value+M)
      return(mu)
    }
  }

  if(fromto==2){
    if(type==1){
      F<--log(1-value)
      return(F)
    }
    if(type==2){
      if(length(value)==1){
      ff<-function(x,value,M){
        d<-(x*(1-exp(-x-M)))/(x+M)
        return(abs(value-d))
      }
      F<-optimize(ff,c(0,100),tol=0.0000001,value=value,M=M)[[1]]
      return(F)
      }
      if(length(value)>1){
        parms<--log(1-value)
        fff<-function(x){
          y<-(x*(1-exp(-x-M)))/(x+M)
          return(sum((value-y)^2))
        }
        Fp<-optim(parms,fff,method="BFGS")
        return(Fp$par)
      }
    }
  }
}

