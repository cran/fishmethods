catch.select<-function(len=NULL, lenmin=NULL, binsize=NULL, peakplus = 1, Linf=NULL,K=NULL,t0=NULL,subobs=FALSE){
if(is.null(lenmin)) stop("lenmin is missing.")
if(is.null(Linf)||is.null(K)||is.null(t0)) stop("Linf, K or t0 is missing.")
  if(binsize==0) stop("binsize must be >0")
  maxlen<-max(len)
  lens<-seq(lenmin,maxlen,binsize)
  datar<-data.frame(lint=c(lens),uint=0,midlen=0,t=0,dt=0,catch=0,logCdt=0,stobs=0,loginvs=0,stest=0)
datar$uint<-datar$lint+binsize
datar$catch[1:c(length(datar[,1]))]<-as.data.frame(table(cut(len, breaks=c(datar$lint,datar$lint[length(datar[,1])]+binsize),include.lowest=FALSE,right=FALSE)))[,2]
if(max(datar$uint)>Linf) warning("Some length intervals are greater than Linf.  These will not be used in the calculations.")
index1<-which(datar$uint<=Linf)
datar$t[index1]<-(((log((datar$lint[index1]-Linf)/-Linf)/-K)+t0)+((log((datar$uint[index1]-Linf)/-Linf)/-K)+t0))/2
datar$dt[index1]<-(1/K)*log((Linf-datar$lint[index1])/(Linf-datar$uint[index1]))
datar$logCdt<-ifelse(datar$catch<=0,0,log(datar$catch/datar$dt))
datar$logCdt<-ifelse(is.infinite(datar$logCdt),0,datar$logCdt)
datar$midlen<-(datar$lint+datar$uint)/2
tempdat<-datar[datar$logCdt>0,]
## Catch Curve
index1<-which(tempdat$logCdt==max(tempdat$logCdt,na.rm=TRUE))+peakplus
outs<-lm(tempdat$logCdt[index1:length(tempdat[,1])]~tempdat$t[index1:length(tempdat[,1])])
p1<-coef(outs)[1]+coef(outs)[2]*tempdat[1:c(index1-1),"t"]
p2<-tempdat[1:c(index1-1),"dt"]*exp(p1)
p3<-tempdat[1:c(index1-1),"catch"]/p2
tempdat$stobs[1:c(index1-1)]<-p3
tempdat$stobs<-ifelse(tempdat$stobs>1,1,tempdat$stobs)
if(subobs==TRUE) tempdat$stobs<-ifelse(tempdat$stobs==1,0.9999,tempdat$stobs)
tempdat$loginvs[1:c(index1-1)]<-ifelse(tempdat$stobs[1:c(index1-1)]<1,log((1/tempdat$stobs[1:c(index1-1)])-1),0)
if(length(tempdat[tempdat[,8]>0 & tempdat[,8]<1,1])==2) warning("Only two observations are available to fit the ogive.")
if(length(tempdat[tempdat[,8]>0 & tempdat[,8]<1,1])==1) stop("There is only one observation available to fit the ogive.")
tempdat$predlogCdt[index1:length(tempdat[,1])]<-coef(outs)[1]+coef(outs)[2]*tempdat[index1:length(tempdat[,1]),"t"]
dodo1<-tempdat[tempdat$stobs>0 & tempdat$stobs<1,]
## Selectivity
outs2<-lm(dodo1$loginvs~dodo1$t)
datar$stest<-1/(1+exp(c(coef(outs2)[1])+c(coef(outs2)[2])*datar$t))
datar$stest<-ifelse(datar$stest<0,0,datar$stest)
datar<-merge(datar,tempdat[c(1,8,11)],by.x="lint",by.y="lint",all.x=TRUE)
datar$logCdt<-ifelse(datar$catch<=0,NA,datar$logCdt)
datar$logCdt<-ifelse(datar$t<=0,NA,datar$logCdt)
datar$stest<-ifelse(datar$t<=0,NA,datar$stest)
estimates<-datar[,c(1,2,3,4,6,7,12,11,10)]
names(estimates)<-c("lower_length","upper_length","mid_length","age_at_midlen","catch","logCdt","predlogCdt","observed_selectivity","estimated_selectivity")
ccparms<-summary(outs)$coefficients[1:2,1:2]
dimnames(ccparms)<-list(c("Intercept","Z"),c("Estimate","Std. Error"))
sparms<-summary(outs2)$coefficients[1:2,1:2]
dimnames(sparms)<-list(c("T1","T2"),c("Estimate","Std. Error"))
ans<-list(estimates=estimates,catchcurve=ccparms,select_ogive=sparms)
return(ans)
}

