vbmultifit<-function(len=NULL,age=NULL,group=NULL,fixed=c(1,1,1),error=1,
        select=1,Linf=c(NULL),K=c(NULL),t0=c(NULL),graph=0,
                control=list(maxiter=10000,minFactor=1/1024,tol=1e-5)){
   if(is.null(len)) 
         stop ("len is missing") 
   if(is.null(age)) 
         stop ("age is missing") 
   if(is.null(group)) 
         stop ("group is missing.") 
   if(length(age)!=length(len)) stop ("Vectors of different lengths")
   if(length(len)!=length(group)) stop ("Vectors of different lengths")
   if(length(age)!=length(group)) stop ("Vectors of different lengths")
  ngroups<-nlevels(as.factor(group)) 
  if(nlevels(as.factor(group))<2) stop("Only two or more groups are allowed.")
  if(select==2 & (is.null(Linf)|is.null(K)|is.null(t0))) stop("User-specified values of Linf, K, and t0 are required")
  cat<-as.data.frame(model.matrix(lm(age~as.factor(group))))
  names(cat)<-levels(group)
x2<-NULL;wgt<-NULL   
x<-as.data.frame(cbind(len,age,cat))
  index<-as.numeric(which(is.na(x),arr.ind=TRUE))[1]
  if(!is.na(index)) x<-x[-index,]
  storeLinf<-NULL; storeK<-NULL;storet0<-NULL
  subcats<-sapply(2:ncol(x),function(y) list(x[,y]))
  subcats<-subcats[c(2:length(subcats),1)]
  subcats[[ngroups+1]]<-trunc(subcats[[ngroups+1]])
  g1<-aggregate(x[,1],subcats,mean)
  names(g1)[ngroups+1]<-"age"
  for(t in 2:ngroups){
    m1<-g1[g1[,t]==1,]
    m1$x2<-NA
    m1$x2[1:c(length(m1$x)-1)]<-m1$x[2:c(length(m1$x))]   
    out1<-lm(x~x2,data=m1,subset=(!is.na(x2)))
    KK<-abs(log(coef(out1)[2]))
    LI<--coef(out1)[1]/(coef(out1)[2]-1)
    dx1<-as.data.frame(cbind(LI-m1$x,m1$age));dx1<-dx1[dx1[,1]>0,]
    t0d<-(coef(lm(log(dx1[,1])~dx1[,2]))[1]-log(LI))/KK
    storeLinf<-c(storeLinf,as.numeric(LI))
    storeK<-c(storeK,as.numeric(KK))
    storet0<-c(storet0,as.numeric(t0d))
  }
  if(select==1){ 
    Lparms<-c(max(storeLinf),-(max(storeLinf)-storeLinf))
    Kparms<-c(max(storeK),max(storeK)-storeK)
    indexs<-which(abs(storet0)==max(abs(storet0)))
    t0parms<-c(storet0[indexs],storet0[indexs]-storet0)
  }
  if(select==2){
    if(length(Linf)!=ngroups) stop("Number of Linfs does not match number of groups")
    if(length(K)!=ngroups) stop("Number of Ks does not match number of groups")
    if(length(t0)!=ngroups) stop("Number of t0s does not match number of groups")
    Lparms<-c(Linf[1],Linf[1]-Linf[2:ngroups])
    Kparms<-c(K[1],K[1]-K[2:ngroups])
    t0parms<-c(t0[1],t0[1]-t0[2:ngroups])
    
  }
  if(error==1) x$wgt<-1
  if(error==2){
    subcats<-sapply(2:ncol(x),function(y) list(x[,y]))
    subcats<-subcats[c(2:length(subcats),1)]
    x<-aggregate(x$len,subcats,mean)
    names(x)<-c(levels(group),"age","len")
    x<-x[,c("len","age",levels(group))]
    x$wgt<-1
    
  }
  if(error==3){
    subcats<-sapply(2:ncol(x),function(y) list(x[,y]))
    subcats<-subcats[c(2:length(subcats),1)]
    d1<-aggregate(x$len,subcats,mean,na.rm=TRUE)
    names(d1)<-c(levels(group),"age","len")
    d2<-aggregate(x$len,subcats,var,na.rm=TRUE)
    names(d2)<-c(levels(group),"age","s2")
    d3<-aggregate(x$len,subcats,length)
    names(d3)<-c(levels(group),"age","n")
    d4<-merge(d1,d2,by.x=c(levels(group),"age"),by.y=c(levels(group),"age"))
    d4<-merge(d3,d4,by.x=c(levels(group),"age"),by.y=c(levels(group),"age"))
    d4$wgt<-d4$n/d4$s2 
    x<-d4[,c("len","age",levels(group),"wgt")]
    if(any(is.na(x$wgt))) stop("At least one age has a single length observation. Need at least two observations to calculate variance." )
  }
  Linfs<-NULL;Ks<-NULL;t0s<-NULL     
  for(pp in 1:length(Lparms)){
    if(pp<length(Lparms)){
      Linfs<-c(Linfs,paste("Linf",pp,"=",round(Lparms[pp],3),",",sep="")) 
      Ks<-c(Ks,paste("K",pp,"=",round(Kparms[pp],3),",",sep=""))
      t0s<-c(t0s,paste("t0",pp,"=",round(t0parms[pp],3),",",sep=""))
    }
    if(pp==length(Lparms)){
      Linfs<-c(Linfs,paste("Linf",pp,"=",round(Lparms[pp],3),sep="")) 
      Ks<-c(Ks,paste("K",pp,"=",round(Kparms[pp],3),sep=""))
      t0s<-c(t0s,paste("t0",pp,"=",round(t0parms[pp],3),sep=""))
    }
  }
  Linfs<-paste(Linfs,collapse="")
  Ks<-paste(Ks,collapse="")
  t0s<-paste(t0s,collapse="")
    #Start here need to change code for fixed parameters
  if(fixed[1]==2){
    pL<-paste("Linf",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
    pL<-paste(c(paste(pL[1:c(ngroups-1)],"+",sep=""),pL[ngroups]),collapse="")
  }
  if(fixed[2]==2){
    pK<-paste("K",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
    pK<-paste(c(paste(pK[1:c(ngroups-1)],"+",sep=""),pK[ngroups]),collapse="")
  }
  if(fixed[3]==2){
    pt0<-paste("t0",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
    pt0<-paste(c(paste(pt0[1:c(ngroups-1)],"+",sep=""),pt0[ngroups]),collapse="")
  }
  if(fixed[1]==1){
    Linfs<-substr(Linfs,1,c(gregexpr(pattern =',',Linfs)[[1]][1]-1))
    pL<-paste("Linf",1,"*",names(x)[3],sep="")
  }
  if(fixed[2]==1){
    Ks<-substr(Ks,1,c(gregexpr(pattern =',',Ks)[[1]][1]-1))
    pK<-paste("K",1,"*",names(x)[3],sep="")
  }
  if(fixed[3]==1){
    t0s<-substr(t0s,1,c(gregexpr(pattern =',',t0s)[[1]][1]-1))
    pt0<-paste("t0",1,"*",names(x)[3],sep="")
  }
  # Full model
  starts<-paste("list(",Linfs,",",Ks,",",t0s,")",sep="",collapse="")
  equat<-paste("len~(",pL,")*","(1-exp(-(",pK,")*(age-(",pt0,"))))",sep="",collapse="")
  Ho<-try(nls(eval(parse(text=equat)),data=x,       
              weights=wgt,start=eval(parse(text=starts)),
              control=control),silent=TRUE)
  if(class(Ho)=="try-error") stop(paste("Ho: ",attributes(Ho)[2],sep=""))
  resid0<-residuals(Ho)
      nlsout<-list(summary(Ho),AIC(Ho),resid0)
      names(nlsout)<-c("results","AIC","residuals")
   
      # Plots
      if(graph==1){
        primcolors<-c("black","red","green3","salmon","blue","purple","orange","gray36","pink")
        par(mfrow=c(1,2))
            cols<-primcolors[1:ngroups]
            getpred<-paste("x$pred<-predict(Ho)",sep="",collapse="")
            eval(parse(text=getpred))
            for(i in 1:ngroups){
              if(i==1){
                temp<-data.matrix(x[,4:c(4+ngroups-2)])
                keep <- apply(temp, 1, function (y){
                  all(y==0)})
                temp1<-x[keep, ]
                temp1<-temp1[order(temp1$age),]
                plot(len~age,data=temp1,lwd=1.5,xlab="Age",ylab="Length",ylim=c(0,max(x$len)),
                     xlim=c(0,max(x$age)),pch=16)
                lines(temp1$pred~temp1$age,col=cols[i])
              }
              if(i>1){
                temp1<-x[x[,3+i-1]==1,]
                temp1<-temp1[order(temp1$age),]
                points(temp1$len~temp1$age,col=cols[i],pch=16)
                lines(temp1$pred~temp1$age,col=cols[i],lwd=1.5)
              }
            }
           legend("topleft", levels(group),bty="n", xpd = TRUE, title="Group",ncol=2, inset = c(0.1,0),x.intersp=0.65,
                 pch=16,cex = 1.0,col=primcolors[1:ngroups])
            getpred<-paste("x$pred<-resid(Ho)",sep="",collapse="")
            eval(parse(text=getpred))
            maxres<-max(abs(x$pred))
            for(i in 1:ngroups){
              if(i==1){
                temp<-data.matrix(x[,4:c(4+ngroups-2)])
                keep <- apply(temp, 1, function (y){
                  all(y==0)})
                temp1<-x[keep, ]
                plot(pred~age,data=temp1,xlab="Age",ylab="Residuals",ylim=c(-maxres,maxres),
                     xlim=c(0,max(x$age)),pch=16)
                abline(h=0)
              }
              if(i>1){
                temp1<-x[x[,3+i-1]==1,]
                points(temp1$pred~temp1$age,col=cols[i],pch=16)
              }
            }
          legend("topleft", levels(group),bty="n", xpd = TRUE, title="Group",ncol=2, inset = c(0.1,0),x.intersp=0.65,
                 pch=16,cex = 1.0,col=primcolors[1:ngroups])
      }#graph==1
  return(nlsout)
}


