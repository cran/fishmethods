mature<-function(cap_age=NULL,mature_age=NULL,age_all_immature=NULL,
age_all_mature=NULL,initial=NULL,nrandoms=1000){
if(length(cap_age)!=length(mature_age)) stop("cap_age and mature_age vector lengths are different")
nparms<-c(age_all_mature-age_all_immature-2)
if(is.null(initial)) stop("At least one starting value in required in initial")
if(any(initial>1)) stop("initial values must be <=1")
if(any(initial<=0)) stop("initial values must be >0")
if(nparms==length(initial)) parms<-initial
if(nparms!=length(initial)){
 if(nparms>length(initial)) parms<-c(initial,rep(initial[length(initial)],nparms-length(initial)))
 if(nparms<length(initial)) parms<-c(initial[1:nparms])
}
if(sum(parms)>1) stop("sum of initial values >1")
parms<-log(parms/(1-sum(parms)))
maxage<-max(max(cap_age),max(mature_age))
minage<-min(min(cap_age),min(mature_age))
nages<-c(maxage-minage+1)
datar<-matrix(0,nrow=maxage,ncol=maxage)
tempdat<-data.frame(acap=cap_age,mage=mature_age)

for(rr in 1:maxage){
  for(cc in 1:maxage){
    tempo<-tempdat[tempdat$acap==rr & tempdat$mage==cc,]
     if(length(tempo[,1])>0) datar[rr,cc]<-length(tempo[,1])
  }
}
colnames(datar)<-c(1:maxage)
rownames(datar)<-c(1:maxage)

datar<-datar[c(c(age_all_immature+1):maxage),c(c(age_all_immature+1):c(age_all_mature-1))]
cutage<-ncol(datar)

solvefor<-function(p){
pp<-exp(p)/(1+sum(exp(p)))
pend<-1/(1+sum(exp(p)))
nll<-matrix(0,nrow=length(datar[,1]),ncol=nparms+1)
for(rr in length(datar[,1]):1){
 if(rr>=cutage){
  for(ll in 1:c(nparms+1)){
    if(ll<cutage) nll[rr,ll]<-newdata[rr,ll]*log(pp[ll])
    if(ll==cutage){
         if(pend<=0) nll[rr,ll]<-newdata[rr,ll]*0
         if(pend>0) nll[rr,ll]<-newdata[rr,ll]*log(pend)
      }  
   }
  }#if rr>=cutage
 if(rr<cutage){
   for(ll in 1:c(nparms+1)){  
      tprob<-1-pend
      if(c(cutage-rr)>1){
        addprob<-sum(pp[c(length(pp)-c(cutage-rr)+2):length(pp)])    
        tprob<-tprob-addprob      
      }   
      if(ll<=rr){
        if(tprob<=0) nll[rr,ll]<-newdata[rr,ll]*0
        if(tprob>0) nll[rr,ll]<-newdata[rr,ll]*log(pp[ll]/tprob)
      }
   }
  } 
}#rr loop
return(-sum(nll,na.rm=TRUE))
}

storep<-matrix(0,nrow=c(nparms+1),ncol=nrandoms)
for(ff in 1:nrandoms){
newdata<-apply(datar,1,function(x){
if(all(x<=0)) rep(0,length(x)) else rmultinom(1,sum(x),x/sum(x))})
newdata<-t(newdata)
outs1<-optim(parms,solvefor,method="Nelder-Mead")
storep[,ff]<-c(exp(outs1$par)/(1+sum(exp(outs1$par))),1/(1+sum(exp(outs1$par))))
}
prms<-apply(storep,1,function(x){mean(x)})
ses<-apply(storep,1,function(x){sd(x)})
#Get expected values
pp<-prms[1:nparms]
pend<-prms[nparms+1]
expected<-matrix(0,nrow=nrow(datar),ncol=nparms+1)
for(rr in length(datar[,1]):1){
  if(rr>=cutage){
    for(ll in 1:c(nparms+1)){
      if(ll<cutage){
        expected[rr,ll]<-sum(datar[rr,])*pp[ll]
      }  
     if(ll==cutage){
        expected[rr,ll]<-sum(datar[rr,])*pend
      }  
    }
  }#if rr>=cutage
  if(rr<cutage){
    for(ll in 1:c(nparms+1)){  
      tprob<-1-pend
      if(c(cutage-rr)>1){
        addprob<-sum(pp[c(length(pp)-c(cutage-rr)+2):length(pp)])    
        tprob<-tprob-addprob      
      } 
      if(ll<=rr){
        if(tprob<=0) expected[rr,ll]<-sum(datar[rr,])*0
        if(tprob>0) expected[rr,ll]<-sum(datar[rr,])*(pp[ll]/tprob)
      }
    }
  }
}#rr loop

colnames(datar)<-c(c(age_all_immature+1):c(age_all_mature-1))
colnames(expected)<-c(c(age_all_immature+1):c(age_all_mature-1))
rownames(expected)<-c(c(age_all_immature+1):maxage)
ests<-data.frame(ML=round(prms,3),SE=round(ses,3))
rownames(ests)<-c(c(age_all_immature+1):c(age_all_mature-1))
ans<-list(estimates=ests,data=datar,expected=expected)
return(ans)
}
