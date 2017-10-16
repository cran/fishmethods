vblrt<-function(len=NULL,age=NULL,group=NULL,error=1,select=1,Linf=c(NULL),
K=c(NULL),t0=c(NULL),plottype=0,
                control=list(maxiter=10000,minFactor=1/1024,tol=1e-5)){
   if(is.null(len)) 
         stop ("len is missing") 
   if(is.null(age)) 
         stop ("age is missing") 
   if(is.null(group)) 
         stop ("group is missing.") 
   group<-as.factor(group)
   if(length(age)!=length(len)) stop ("Vectors of different lengths")
    ngroups<-nlevels(group) 
    if(nlevels(group)<2) stop("Only two or more groups are allowed.")
   if(select==2 & (is.null(Linf)|is.null(K)|is.null(t0))) stop("User-specified values of Linf, K, and t0 are required")
    cat<-as.data.frame(model.matrix(lm(age~as.factor(group))))
    x2<-NULL;wgt<-NULL  
     names(cat)<-levels(group)
     	x<-as.data.frame(cbind(len,age,cat))
      index<-as.numeric(which(is.na(x),arr.ind=TRUE))[1]
     	if(!is.na(index)) x<-x[-index,]
         
      if(select==1){ 
           storeLinf<-NULL; storeK<-NULL;storet0<-NULL
           subcats<-sapply(2:ncol(x),function(y) list(x[,y]))
           subcats<-subcats[c(2:length(subcats),1)]
           subcats[[ngroups+1]]<-trunc(subcats[[ngroups+1]])
           g1<-aggregate(x[,1],subcats,mean)
           names(g1)[ngroups+1]<-"age"

           for(t in 2:ngroups){
             m1<-g1[g1[,t]==1,]
             temp<-m1[,c(ncol(m1)-1,ncol(m1))];names(temp)<-c("age","len")
             m1$age<-m1$age+1
             m1<-merge(m1,temp,by.x="age",by.y="age",all.x=TRUE,all.y=TRUE)
             m1<-m1[-unique(which(is.na(m1),arr.ind=TRUE)[,1]),]   
             out1<-lm(len~x,data=m1)
             KK<-as.numeric(abs(log(coef(out1)[2])))
             LI<-coef(out1)[1]/(1-coef(out1)[2])
             dx1<-as.data.frame(cbind(LI-m1$x,m1$age));dx1<-dx1[dx1[,1]>0,]
             t0d<-(coef(lm(log(dx1[,1])~dx1[,2]))[1]-log(LI))/KK
             storeLinf<-c(storeLinf,as.numeric(LI))
             storeK<-c(storeK,as.numeric(KK))
             storet0<-c(storet0,as.numeric(t0d))
           } 
           qL<-quantile(storeLinf,prob=c(0.1,0.90))
           storeLinf<-ifelse(storeLinf<qL[1]|storeLinf>qL[2],median(storeLinf),storeLinf)
           qK<-quantile(storeK,prob=c(0.1,0.90))
           storeK<-ifelse(storeK<qK[1]|storeK>qK[2],median(storeK),storeK)
           qt<-quantile(storet0,prob=c(0.1,0.90))
           storet0<-ifelse(storet0<qt[1]|storet0>qt[2],median(storet0),storet0)
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
      pL<-paste("Linf",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
      pL<-paste(c(paste(pL[1:c(ngroups-1)],"+",sep=""),pL[ngroups]),collapse="")
      pK<-paste("K",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
      pK<-paste(c(paste(pK[1:c(ngroups-1)],"+",sep=""),pK[ngroups]),collapse="")
      pt0<-paste("t0",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
      pt0<-paste(c(paste(pt0[1:c(ngroups-1)],"+",sep=""),pt0[ngroups]),collapse="")
      
      
    # Full model
      starts<-paste("list(",Linfs,",",Ks,",",t0s,")",sep="",collapse="")
      equat<-paste("len~(",pL,")*","(1-exp(-(",pK,")*(age-(",pt0,"))))",sep="",collapse="")
      Ho<-try(nls(eval(parse(text=equat)),data=x,       
        	 weights=wgt,start=eval(parse(text=starts)),
             control=control),silent=TRUE)
      if(class(Ho)=="try-error") stop(paste("Ho: ",attributes(Ho)[2],sep=""))
      resid0<-residuals(Ho)
      AICo<-AIC(Ho)
    # H1
      Linf1<-substr(Linfs,1,c(gregexpr(pattern =',',Linfs)[[1]][1]-1))
      pK<-paste("K",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
      pK<-paste(c(paste(pK[1:c(ngroups-1)],"+",sep=""),pK[ngroups]),collapse="")
      pt0<-paste("t0",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
      pt0<-paste(c(paste(pt0[1:c(ngroups-1)],"+",sep=""),pt0[ngroups]),collapse="")
      pL<-paste("Linf",1,"*",names(x)[3],sep="")
      starts<-paste("list(",Linf1,",",Ks,",",t0s,")",sep="",collapse="")
      equat<-paste("len~(",pL,")*","(1-exp(-(",pK,")*(age-(",pt0,"))))",sep="",collapse="")
      H1<-try(nls(eval(parse(text=equat)),data=x,        
		       weights=wgt,start=eval(parse(text=starts)),control=control),silent=TRUE)
      if(class(H1)=="try-error") stop(paste("H1: ",attributes(H1)[2],sep=""))
	    resid1<-residuals(H1)
	    AIC1<-AIC(H1)
	#H2
	    K2<-substr(Ks,1,c(gregexpr(pattern =',',Ks)[[1]][1]-1))
	    pK<-paste("K",1,"*",names(x)[3],sep="")
	    pL<-paste("Linf",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
	    pL<-paste(c(paste(pL[1:c(ngroups-1)],"+",sep=""),pL[ngroups]),collapse="")
	    pt0<-paste("t0",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
	    pt0<-paste(c(paste(pt0[1:c(ngroups-1)],"+",sep=""),pt0[ngroups]),collapse="")
	    
	    starts<-paste("list(",Linfs,",",K2,",",t0s,")",sep="",collapse="")
	    equat<-paste("len~(",pL,")*","(1-exp(-(",pK,")*(age-(",pt0,"))))",sep="",collapse="")
 	    H2<-try(nls(eval(parse(text=equat)),data=x,       
 	         weights=wgt,start=eval(parse(text=starts)),
 	         control=control),silent=TRUE)
 	    if(class(H2)=="try-error") stop(paste("H2: ",attributes(H2)[2],sep=""))
     resid2<-residuals(H2)
     AIC2<-AIC(H2)
               
# H3      
     pL<-paste("Linf",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
     pL<-paste(c(paste(pL[1:c(ngroups-1)],"+",sep=""),pL[ngroups]),collapse="")
     pK<-paste("K",1:ngroups,"*",names(x)[3:c(ncol(x)-1)],sep="")
     pK<-paste(c(paste(pK[1:c(ngroups-1)],"+",sep=""),pK[ngroups]),collapse="")
     t01<-substr(t0s,1,c(gregexpr(pattern =',',t0s)[[1]][1]-1))
     pt0<-paste("t0",1,"*",names(x)[3],sep="")
     starts<-paste("list(",Linfs,",",Ks,",",t01,")",sep="",collapse="")
     equat<-paste("len~(",pL,")*","(1-exp(-(",pK,")*(age-(",pt0,"))))",sep="",collapse="")
     H3<-try(nls(eval(parse(text=equat)),data=x,       
    	         weights=wgt,start=eval(parse(text=starts)),
    	         control=control),silent=TRUE)
     if(class(H3)=="try-error") stop(paste("H3: ",attributes(H3)[2],sep=""))
     resid3<-residuals(H3)
     AIC3<-AIC(H3)
     #H4
     t01<-substr(t0s,1,c(gregexpr(pattern =',',t0s)[[1]][1]-1))
     Linf1<-substr(Linfs,1,c(gregexpr(pattern =',',Linfs)[[1]][1]-1))
     K1<-substr(Ks,1,c(gregexpr(pattern =',',Ks)[[1]][1]-1))
     pL<-paste("Linf",1,"*",names(x)[3],sep="")
     pK<-paste("K",1,"*",names(x)[3],sep="")
     pt0<-paste("t0",1,"*",names(x)[3],sep="")
     starts<-paste("list(",Linf1,",",K1,",",t01,")",sep="",collapse="")
        equat<-paste("len~(",pL,")*","(1-exp(-(",pK,")*(age-(",pt0,"))))",sep="",collapse="")
    H4<-try(nls(eval(parse(text=equat)),data=x,       
                weights=wgt,start=eval(parse(text=starts)),
                control=control),silent=TRUE)
    if(class(H4)=="try-error") stop(paste("H4: ",attributes(H4)[2],sep=""))
         resid4<-residuals(H4)
         AIC4<-AIC(H4)
  	 RSS<-c(sum(residuals(Ho)^2),sum(residuals(H1)^2),sum(residuals(H2)^2),
             sum(residuals(H3)^2),sum(residuals(H4)^2))
 	   N<-length(residuals(Ho))
  	 X<-round(c(-N*log(RSS[1]/RSS[2]),-N*log(RSS[1]/RSS[3]),-N*log(RSS[1]/RSS[4]),
             -N*log(RSS[1]/RSS[5])),2)
  	 df<-c(length(coef(Ho))-length(coef(H1)),length(coef(Ho))-length(coef(H2)),
       	length(coef(Ho))-length(coef(H3)),length(coef(Ho))-length(coef(H4)))
  	 p<-round(1-pchisq(X,df),3)
     AICC<-c(round(AICo,2),round(AIC1,2),round(AIC2,2),round(AIC3,2),round(AIC4,3))
      labs<-c("Ho","H1","H2","H3","H4")
      Llabs<-paste(c(paste("Linf",1:c(ngroups-1),"=",sep=""),paste("Linf",ngroups,sep="")),collapse="")
      Klabs<-paste(c(paste("K",1:c(ngroups-1),"=",sep=""),paste("K",ngroups,sep="")),collapse="")
      t0labs<-paste(c(paste("t0",1:c(ngroups-1),"=",sep=""),paste("t0",ngroups,sep="")),collapse="")
      
      hyp<-c(Llabs,Klabs,t0labs,paste(Llabs,",",Klabs,",",t0labs,sep=""))
      labels<-c("Ho vs H1","Ho vs H2","Ho vs H3","Ho vs H4")
      compout<-data.frame(tests=labels,hypothesis=hyp,chisq=X,df=df,p=p)
      rss<-as.data.frame(cbind(labs,RSS,AICC));names(rss)<-c("model","rss","AIC")
      residuals_all<-as.data.frame(cbind(resid0,resid1,resid2,resid3,resid4))
      nlsout<-list(compout,summary(Ho),summary(H1),summary(H2),summary(H3), summary(H4),
             rss,residuals_all)
      names(nlsout)<-c("results",c(paste("model",labs)),"rss","residuals")
    # Plot observed versus predicted
    if (plottype>0){
      if(plottype==1){
 	      	par(mfrow=c(3,2))
          primcolors<-c("black","red","green3","salmon","blue","purple","orange","gray36","pink")
          cols<-primcolors[1:ngroups]
         for(plotter in 1:5){
            if(plotter==1) lalab<-"Ho"
            if(plotter>1) lalab<-paste("H",plotter-1,sep="",collapse="")
            getpred<-paste("x$pred<-predict(",lalab,")",sep="",collapse="")
            eval(parse(text=getpred))
		       for(i in 1:ngroups){
		        if(i==1){
		        temp<-data.matrix(x[,4:c(4+ngroups-2)])
		        keep <- apply(temp, 1, function (y){
		        all(y==0)})
		        temp1<-x[keep, ]
		        temp1<-temp1[order(temp1$age),]
		        plot(len~age,data=temp1,main=c(paste(lalab," Model",sep="")),lwd=1.5,xlab="Age",ylab="Length",ylim=c(0,max(x$len)),
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
      
                        
     } #Plotter loop
    plot(0,0, type = "n", bty = "n", xaxt = "n", yaxt = "n",col.lab="transparent")
    legend("topleft", levels(group),bty="n", xpd = TRUE, title="Group",ncol=2, inset = c(0.1,0),x.intersp=0.65,
    pch=16,cex = 1.5,col=primcolors[1:ngroups])
      par(mfrow=c(1,1))
		}#plottype==1
   if (plottype==2){
     par(mfrow=c(3,2))
     # Ho model
     primcolors<-c("black","red","green3","salmon","blue","purple","orange","gray36","pink")
     cols<-primcolors[1:ngroups]
  for(plotter in 1:5){
    if(plotter==1) lalab<-"Ho"
    if(plotter>1) lalab<-paste("H",plotter-1,sep="",collapse="")
     getpred<-paste("x$pred<-resid(",lalab,")",sep="",collapse="")
     eval(parse(text=getpred))
     if(plotter==1) maxres<-max(abs(x$pred))
     if(plotter>1){
       maxres<-max(maxres,max(abs(x$pred)))
     }
     for(i in 1:ngroups){
       if(i==1){
         temp<-data.matrix(x[,4:c(4+ngroups-2)])
         keep <- apply(temp, 1, function (y){
           all(y==0)})
         temp1<-x[keep, ]
         plot(pred~age,data=temp1,main=c(paste(lalab," Model",sep="")),xlab="Age",ylab="Residuals",ylim=c(-maxres,maxres),
              xlim=c(0,max(x$age)),pch=16)
         abline(h=0)
       }
       if(i>1){
         temp1<-x[x[,3+i-1]==1,]
         points(temp1$pred~temp1$age,col=cols[i],pch=16)
       }
     }
     }#plotter loop
     plot(0,0, type = "n", bty = "n", xaxt = "n", yaxt = "n",col.lab="transparent")
     legend("topleft", levels(group),bty="n", xpd = TRUE, title="Group",ncol=2, inset = c(0.1,0),x.intersp=0.65,
            pch=16,cex = 1.5,col=primcolors[1:ngroups])
     par(mfrow=c(1,1))
  }#plottype==2
}#plottype>0
  return(nlsout)
}
