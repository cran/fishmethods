fm_telemetry<-function(filetype=c(1), caphistory=NULL,
                Fdesign=NULL, Mdesign=NULL, Pdesign=NULL,
                whichlivecells=NULL,whichdeadcells=NULL,constant=1e-14,
                initial=NULL,invtol=1e-44,
                control=list(reltol=1e-8,maxit=1000000)){
if(filetype==1) datar<-caphistory
if(filetype==2)  datar<-scan(caphistory,what="character")
datar<-sort(datar,decreasing=TRUE)
if(!all(nchar(datar[1])==nchar(datar))) stop("histories are not the same length")
index<-which(datar=="NA",arr.ind=TRUE)
if(length(index)>0) stop(paste("There are NAs in row(s) ", paste0('"', paste(as.character(index[1]), collapse="\", \""), '"'),sep=""))

if(any(is.numeric(Fdesign))) stop("Fdesign must be a character vector")
if(any(is.numeric(Mdesign))) stop("Mdesign must be a character vector")
if(any(is.numeric(Pdesign))) stop("Pdesign must be a character vector")

# Check if occassion are coded incorrectly
nocc<-nchar(datar[1])
occs<-1:nocc
warns<-NULL
   for(i in 1:length(datar)){
     for(occ in 1:nocc){
      temp<-substr(datar[i],occs[occ],occs[occ])
       if(!temp %in% c("0","1","D","E")) warns<-rbind(warns,i)
     }
   }
  if(length(warns)>0) stop(paste("Capture History: There are coding errors in row(s) ", paste0('"', paste(as.character(warns), collapse="\", \""), '"'),sep=""))
tempo<-NULL
for(i in 1:nocc) tempo<-paste(tempo,"0",sep="") 
datar<-datar[datar!=tempo]

#################### Setup Fdesign ####################################
Design<-matrix(0,nrow=c(nocc-1),4)
Design[,1]<-1:c(nocc-1)
pardef<-strsplit(Fdesign,"*",fixed=TRUE)
for(t in 1:length(pardef)){
  tempo<-pardef[[t]]
  #single or sequence of parameters
  if(length(tempo)==1){
    eval(parse(text=paste("seqs<-",tempo[1],sep="")))
    index<-which(Design[,1] %in% c(seqs))
    first<-c(max(Design[,2])+1)
    remain<-length(seqs)-1
    if(remain>0) dtem<-c(first,c(first+c(1:c(length(seqs)-1)))) else dtem<-first
    Design[index,2]<-dtem
  }
  #Begin cells subsitution
  if(length(tempo)>1){
    eval(parse(text=paste("seqs<-",tempo[1],sep="")))
    index<-which(Design[,1] %in% c(seqs))
    first<-c(max(Design[,2])+1)
    remain<-length(seqs)-1
    if(remain>0) dtem<-c(first,c(first+c(1:c(length(seqs)-1)))) else dtem<-first
    Design[index,2]<-dtem
    eval(parse(text=paste("subcells<-",tempo[2],sep="")))
    index<-which(Design[,1] %in% c(subcells))
    if(length(subcells)>length(dtem)) dtem<-c(dtem,rep(dtem[length(dtem)],c(length(subcells)-length(dtem))))
    if(length(subcells)<length(dtem)) dtem<-dtem[1:length(subcells)]
    Design[index,2]<-dtem
  }   
}    
if(Design[1,2]==0) stop("Occasion 1 is not included in the F design")
if(Design[1,2]>0) for(t in 1:length(Design[,1])){if(Design[t,2]==0) Design[t,2]<-Design[t-1,2]}
maxparm<-max(Design[,2])

################## Setup M design ########################
pardef<-strsplit(Mdesign,"*",fixed=TRUE)
for(t in 1:length(pardef)){
  tempo<-pardef[[t]]
  #single or sequence of parameters
  if(length(tempo)==1){
    eval(parse(text=paste("seqs<-",tempo[1],sep="")))
    index<-which(Design[,1] %in% c(seqs))
    first<-max(maxparm+1,max(Design[,3])+1)
    remain<-length(seqs)-1
    if(remain>0) dtem<-c(first,c(first+c(1:c(length(seqs)-1)))) else dtem<-first
    Design[index,3]<-dtem
  }
  
  #Begin cells subsitution
  if(length(tempo)>1){
    eval(parse(text=paste("seqs<-",tempo[1],sep="")))
    index<-which(Design[,1] %in% c(seqs))
    first<-max(maxparm+1,max(Design[,3])+1)
    remain<-length(seqs)-1
    if(remain>0) dtem<-c(first,c(first+c(1:c(length(seqs)-1)))) else dtem<-first
    Design[index,3]<-dtem
    eval(parse(text=paste("subcells<-",tempo[2],sep="")))
    index<-which(Design[,1] %in% c(subcells))
    if(length(subcells)>length(dtem)) dtem<-c(dtem,rep(dtem[length(dtem)],c(length(subcells)-length(dtem))))
    if(length(subcells)<length(dtem)) dtem<-dtem[1:length(subcells)]
    Design[index,3]<-dtem
  }   
}    
if(Design[1,3]==0) stop("Occasion 1 is not included in the M design")
if(Design[1,3]>0) for(t in 1:length(Design[,1])){if(Design[t,3]==0) Design[t,3]<-Design[t-1,3]}
maxparm<-max(Design[,3])

########################### Setup Prob Detection ##############################################

pardef<-strsplit(Pdesign,"*",fixed=TRUE)
for(t in 1:length(pardef)){
  tempo<-pardef[[t]]
  #single or sequence of parameters
  if(length(tempo)==1){
    eval(parse(text=paste("seqs<-",tempo[1],sep="")))
    seqs<-seqs-1
    index<-which(Design[,1] %in% c(seqs))
    first<-max(maxparm+1,max(Design[,4])+1)
    remain<-length(seqs)-1
    if(remain>0) dtem<-c(first,c(first+c(1:c(length(seqs)-1)))) else dtem<-first
    Design[index,4]<-dtem
  }
  
  #Begin cells subsitution
  if(length(tempo)>1){
    eval(parse(text=paste("seqs<-",tempo[1],sep="")))
    seqs<-seqs-1
    index<-which(Design[,1] %in% c(seqs))
    first<-max(maxparm+1,max(Design[,4])+1)
    remain<-length(seqs)-1
    if(remain>0) dtem<-c(first,c(first+c(1:c(length(seqs)-1)))) else dtem<-first
    Design[index,4]<-dtem
    eval(parse(text=paste("subcells<-",tempo[2],sep="")))
    subcells<-subcells-1
    index<-which(Design[,1] %in% c(subcells))
    if(length(subcells)>length(dtem)) dtem<-c(dtem,rep(dtem[length(dtem)],c(length(subcells)-length(dtem))))
    if(length(subcells)<length(dtem)) dtem<-dtem[1:length(subcells)]
    Design[index,4]<-dtem
  }   
}    
if(Design[1,4]>0) for(t in 1:length(Design[,1])){if(Design[t,4]==0) Design[t,4]<-Design[t-1,4]}

Fp<-rep(log(exp(-initial[1])/(1-exp(-initial[1]))),length(unique(Design[,2])))
Mp<-rep(log(exp(-initial[2])/(1-exp(-initial[2]))),length(unique(Design[,3]))) 
pp<-rep(log(initial[3]/(1-initial[3])),length(unique(Design[,4])))
parms<-c(Fp,Mp,pp)
################ How many Release Cohorts? ######################
cohorts<-0
Rs<-0
if(any(substr(datar,1,1)%in% c("1"))){
  cohorts<-cohorts+1
  Rs<-sum(substr(datar,1,1)=="1")
  locate<-1
 }

for(k in 2:c(nocc-1)){
  val<-NULL
  for(t in 1:c(k-1)) val<-paste(val,"0",sep="")
  val<-paste(val,c("1"),sep="")
  if(any(substr(datar,1,k) %in% val)){
    cohorts<-cohorts+1
    locate<-c(locate,k)
    Rs<-c(Rs,sum(substr(datar,1,k) %in% val))
  }
}  
############ Generate reduced m-matrix ###################################
live<-matrix(0,nrow=c(nocc-1),ncol=c(nocc))
dead<-matrix(0,nrow=c(nocc-1),ncol=c(nocc))
#Alive (set D,E to zero)
datar1<-datar
 for(i in 1:length(datar)){
   for(j in 1:nocc){
     if(substr(datar1[i],j,j) %in% c("D","E")){
      substr(datar1[i],j,j)<-"0"
     }
   }
 }

 for(ind in 1:length(datar1)){
   fish<-datar1[ind]
   for(i in 1:c(nchar(fish)-1)){
     for(j in c(i+1):c(nchar(fish))){ 
        if(substr(fish,i,i)=="1" & substr(fish,j,j)=="1"){
          live[i,j]<-live[i,j]+1
          break
        }
     }
   }
 }
    #Dead
datar2<-datar
 for(ind in 1:length(datar2)){
    fish<-datar2[ind]
    for(i in 1:c(nchar(fish)-1)){
      for(j in c(i+1):c(nchar(fish))){ 
        if(substr(fish,i,i)=="1" & substr(fish,j,j)=="1") break
        if(substr(fish,i,i)=="1" & substr(fish,j,j)=="D"){
          dead[i,j]<-dead[i,j]+1
          break
         }
        }
      }
 } 
   live[locate,1]<-Rs
   dead[locate,1]<-Rs

   #Censor fish
   remRs<-NULL
   datar2<-datar
   for(ind in 1:length(datar2)){
     fish<-datar2[ind]
     for(i in 1:c(nchar(fish)-1)){
         if(substr(fish,i,i)=="1" & substr(fish,i+1,i+1)=="E"){
           remRs<-c(remRs,i)
           break
         }
       }
     }
    
   setR<-as.data.frame(table(remRs))
   setR[,1]<-as.numeric(as.character(setR[,1]))
   live[c(setR[,1]),1]<-live[c(setR[,1]),1]-setR[,2]
   
   ns<-colSums(live[,1:ncol(live)])
    for(i in 2:nrow(live)){
    live[i,1]<-live[i,1]+ns[i]
    dead[i,1]<-dead[i,1]+ns[i]
    }
   never_seen<-live[,1]-rowSums(live[,2:ncol(live)])-rowSums(dead[,2:ncol(dead)])
   marray<-list(live=live,dead=dead,never_seen=never_seen)
  ######################### Create occasions likelihood cells
  liveschedule<-data.frame(occ=1:c(nocc-1),firstcell=2:nocc,lastcell=nocc)
  if(!is.null(whichlivecells[[1]])){
      cn<-length(whichlivecells)
      for(i in 1:cn){
        ted<-whichlivecells[[i]]
        if(ted[2]<ted[1]) stop("end occasion less than start occasion")
        addm<-ted[3]+c(liveschedule[c(ted[1]:ted[2]),2]-1)
        addm[addm>nocc]<-nocc
        liveschedule[c(ted[1]:ted[2]),3]<-addm
      }
  }
  deadschedule<-data.frame(occ=1:c(nocc-1),firstcell=2:nocc,lastcell=nocc)
  if(!is.null(whichdeadcells[[1]])){
    cn<-length(whichdeadcells)
    for(i in 1:cn){
      ted<-whichdeadcells[[i]]
      if(ted[2]<ted[1]) stop("end occasion less than start occasion")
      addm<-ted[3]+c(deadschedule[c(ted[1]:ted[2]),2]-1)
      addm[addm>nocc]<-nocc
      deadschedule[c(ted[1]:ted[2]),3]<-addm
    }
  }


#Binomial Constant
bincof<-0
#Total N
for(rr in 1:c(nocc-1)){
  nfact_in<-marray[[1]][rr,1]
  nfact_out<-0
  if(nfact_in>=2){
    for(kk in 2:nfact_in) nfact_out<-nfact_out+log(kk)
  }
  bincof<-bincof+nfact_out
}

#Never_seen
for(rr in 1:c(nocc-1)){
  nfact_in<-marray[[3]][rr]
  nfact_out<-0
  if(nfact_in>=2){
    for(kk in 2:nfact_in) nfact_out<-nfact_out+log(kk)
  }
  bincof<-bincof-1.0*nfact_out
}
#live and dead
for(rr in 1:c(nocc-1)){
  for(cc in c(rr+1):nocc){
    nfact_in<-marray[[1]][rr,cc]
    nfact_out<-0
    if(nfact_in>=2){
      for(kk in 2:nfact_in) nfact_out<-nfact_out+log(kk)
    }
    bincof<-bincof-1.0*nfact_out
  }
}
for(rr in 1:c(nocc-1)){
  for(cc in c(rr+1):nocc){
    nfact_in<-marray[[2]][rr,cc]
    nfact_out<-0
    if(nfact_in>=2){
      for(kk in 2:nfact_in) nfact_out<-nfact_out+log(kk)
    }
    bincof<-bincof-1.0*nfact_out
  }
}
################### Fit Model #################################
 fit.obs<-function(x){
  #####calc_F_vector
  Fvector<-rep(0,c(nocc-1))
  Mvector<-rep(0,c(nocc-1))
  Pvector<-rep(0,c(nocc-1))
    Fvector<-1/(1+exp(-x[c(Design[,2])]))
    Mvector<-1/(1+exp(-x[c(Design[,3])]))
    Pvector<-1/(1+exp(-x[c(Design[,4])]))
   
  #### Calculate Survival Matrix #######################
  s<-array(1,dim=c(c(nocc-1),c(nocc)))
  cnt<-0
  for (t in 1:c(nocc-1)){
    for (y in as.numeric(c(cnt+1)):c(nocc)){
      if(t==y){s[t,y]<-1}
      if(t!=y){
        s[t,y]<-Fvector[y-1]*Mvector[y-1]    
      }
    }		   
    cnt<-cnt+1
  }
  cnt<-0
  s_prob<-array(1,dim=c(nocc-1,nocc))
  for (t in 1:c(nocc-1)){
    looper<-0
    for (y in as.numeric(cnt+1):c(nocc)){
      probs<-1
      for(a in as.numeric(y-looper):y){
        probs<-probs*s[t,a]
      }
      s_prob[t,y]<-probs
      looper<-looper+1
    }		   
    cnt<-cnt+1
  }
  ## Probability of detection matrix ##
  pdetect<-c(1,Pvector)
  p_prob<-array(1,dim=c(c(nocc-1),c(nocc)))
  for (t in 1:c(nocc-1)){
    for (y in as.numeric(t+1):c(nocc)){
      if(t==y) p_prob[t,y]<-1
      if(t!=y) p_prob[t,y]<-pdetect[y]
    }		   
  }
  minusp<-array(1,dim=c(c(nocc-1),c(nocc)))
  for (t in 1:c(nocc-2)){
    for (y in as.numeric(t+2):c(nocc)){
      if(t==y) minusp[t,y]<-1
      if(t!=y) minusp[t,y]<-1-pdetect[y-1]
    }		   
  }
  minusp_prob<-array(1,dim=c(nocc-1,nocc))
  for(ee in 1:c(nocc-1)) minusp_prob[ee,]<-cumprod(minusp[ee,])
  prob_alive<-s_prob*p_prob*minusp_prob
 
  # Dead fish
  c_prob<-array(1,dim=c(nocc-1,nocc))
  for (t in 1:c(nocc-1)){
    for (y in as.numeric(t+1):c(nocc)){
      if(t==y) c_prob[t,y]<-1
      if(t!=y) {
        reds<--log(Mvector[y-1])/(-log(Fvector[y-1])+(-log(Mvector[y-1])))
        reds<-ifelse(is.nan(reds),0,reds)
        c_prob[t,y]<-reds
      }
    }		   
  }
  minus_s<-1-s
  s_dead_prob<-array(1,dim=c(nocc-1,nocc))
  for (t in 1:c(nocc-1)){
      s_dead_prob[t,c(t+1):nocc]<-s_prob[t,t:c(nocc-1)]
    }
prob_dead<-c_prob*p_prob*minusp_prob*minus_s*s_dead_prob
######## Calculate Likelihoods ########################
    LLA<-0;LLNS<-0
    #alive
      tempo1<-marray[[1]]
      tempo2<-marray[[2]]
       temp_ll1<-matrix(0,nrow=nrow(tempo1),ncol=nocc)
      temp_ll2<-matrix(0,nrow=nrow(tempo2),ncol=nocc)
      
      for(rr in 1:c(nocc-1)){
        first<-liveschedule[rr,2]
        last<-liveschedule[rr,3]
          for(cc in first:last){
           if(tempo1[rr,cc]>0) temp_ll1[rr,cc]<-tempo1[rr,cc]*log(max(constant,prob_alive[rr,cc]))
         }
      }
      for(rr in 1:c(nocc-1)){
        first<-deadschedule[rr,2]
        last<-deadschedule[rr,3]
        for(cc in first:last){
          if(tempo2[rr,cc]>0) temp_ll2[rr,cc]<-tempo2[rr,cc]*log(max(constant,prob_dead[rr,cc]))
        }
      }
      LLA<-sum(temp_ll1)+sum(temp_ll2)
  # Never Seen
      tempo3<-marray[[3]]
      sumdo<-NULL
      temp_ll3<-rep(0,nocc-1)
      for(rr in 1:c(nocc-1)){
        sumdo<-0
        for(cc in c(rr+1):nocc) sumdo<-sumdo+prob_dead[rr,cc]+prob_alive[rr,cc]
        if(tempo3[rr]>0) temp_ll3[rr]<-tempo3[rr]*log(max(constant,1-sumdo))
      } 
       LLNS<-sum(temp_ll3)
       -(LLA+LLNS+bincof)
       
}# fit.obs  
results<-try(optim(parms, fit.obs,method="BFGS",hessian=TRUE,control=control),silent=TRUE)
if(class(results)!="try-error"){
  Sf<-NULL;Sm<-NULL;P<-NULL
  Sf<-1/(1+exp(-results$par[c(Design[,2])]))
  Sm<-1/(1+exp(-results$par[c(Design[,3])]))
  P<-1/(1+exp(-results$par[c(Design[,4])]))
   
  FF<--log(Sf)
  M<--log(Sm)
  var1<-try(diag(solve(results$hessian,tol=invtol)),silent=TRUE)
   FSE<-NULL
   MSE<-NULL
   PSE<-NULL
   Sfvar<-NULL;Smvar<-NULL;Pvar<-NULL
   if(class(var1)!="try-error"){
   # var1<-ifelse(var1<0,0,var1)
    for(t in 1:c(nocc-1)){
      Sfvar[t]<-(var1[Design[t,2]]*exp(2.0*results$par[Design[t,2]]))/(1+exp(results$par[Design[t,2]]))^4
      Smvar[t]<-(var1[Design[t,3]]*exp(2*results$par[Design[t,3]]))/(1+exp(results$par[Design[t,3]]))^4
      Pvar[t]<-(var1[Design[t,4]]*exp(2*results$par[Design[t,4]]))/(1+exp(results$par[Design[t,4]]))^4
    }
    FSE<-sqrt(Sfvar/Sf^2)
    MSE<-sqrt(Smvar/Sm^2)
    PSE<-sqrt(Pvar)
    }
  if(class(var1)=="try-error"){
    FSE<-rep(NA,length(FF))
    MSE<-rep(NA,length(M))
    PSE<-rep(NA,length(P))
  }

parmtable<-data.frame(Occasion=c(1:nocc),F=0,FSE=0,M=0,MSE=0,p=0,pSE=0)
  parmtable[1:c(nocc-1),2]<-FF;parmtable[1:c(nocc-1),3]<-FSE
  parmtable[1:c(nocc-1),4]<-M;parmtable[1:c(nocc-1),5]<-MSE
  parmtable[2:c(nocc),6]<-P;parmtable[2:c(nocc),7]<-PSE
  
############## Calculate Expected and LIklihoods
#### Calculate Survival Matrix #######################
s<-array(1,dim=c(c(nocc-1),c(nocc)))
cnt<-0
for (t in 1:c(nocc-1)){
  for (y in as.numeric(c(cnt+1)):c(nocc)){
    if(t==y){s[t,y]<-1}
    if(t!=y){
      s[t,y]<-exp(-parmtable[y-1,2]-parmtable[y-1,4])    
    }
  }		   
  cnt<-cnt+1
}
cnt<-0
s_prob<-array(1,dim=c(nocc-1,nocc))
for (t in 1:c(nocc-1)){
  looper<-0
  for (y in as.numeric(cnt+1):c(nocc)){
    probs<-1
    for(a in as.numeric(y-looper):y){
      probs<-probs*s[t,a]
    }
    s_prob[t,y]<-probs
    looper<-looper+1
  }		   
  cnt<-cnt+1
}
## Probability of detection matrix ##
p_prob<-array(1,dim=c(c(nocc-1),c(nocc)))
for (t in 1:c(nocc-1)){
  for (y in as.numeric(t+1):c(nocc)){
    if(t==y) p_prob[t,y]<-1
    if(t!=y) p_prob[t,y]<-parmtable[y,6]
  }		   
}
minusp<-array(1,dim=c(c(nocc-1),c(nocc)))
for (t in 1:c(nocc-2)){
  for (y in as.numeric(t+2):c(nocc)){
    if(t==y) minusp[t,y]<-1
    if(t!=y) minusp[t,y]<-1-parmtable[y-1,6]
  }		   
}
minusp_prob<-array(1,dim=c(nocc-1,nocc))
for(y in 1:c(nocc-1)) minusp_prob[y,]<-cumprod(minusp[y,])
prob_alive<-s_prob*p_prob*minusp_prob
prob_alive[lower.tri(prob_alive, diag = TRUE)]<-0
# Dead fish
c_prob<-array(1,dim=c(nocc-1,nocc))
for (t in 1:c(nocc-1)){
  for (y in as.numeric(t+1):c(nocc)){
    if(t==y) c_prob[t,y]<-1
    if(t!=y) c_prob[t,y]<-ifelse(is.nan(parmtable[y-1,4]/(parmtable[y-1,2]+parmtable[y-1,4])),0,parmtable[y-1,4]/(parmtable[y-1,2]+parmtable[y-1,4]))
  }		   
}
minus_s<-1-s
s_dead_prob<-array(1,dim=c(nocc-1,nocc))
for (t in 1:c(nocc-1)){
  s_dead_prob[t,c(t+1):nocc]<-s_prob[t,t:c(nocc-1)]
}
prob_dead<-c_prob*p_prob*minusp_prob*minus_s*s_dead_prob
prob_dead[lower.tri(prob_dead, diag = TRUE)]<-0

######## Calculate Expected Alive and Dead ########################
tempo1<-marray[[1]]
tempo2<-marray[[2]]
exp_alive<-matrix(0,nrow=nrow(marray[[1]]),ncol=ncol(marray[[1]]))
exp_dead<-matrix(0,nrow=nrow(marray[[1]]),ncol=ncol(marray[[1]]))
exp_ns<-rep(0,nocc-1)
for(rr in 1:c(nocc-1)){
  for(cc in c(rr+1):nocc){
      exp_alive[rr,cc]<-tempo1[rr,1]*prob_alive[rr,cc]
      exp_dead[rr,cc]<-tempo2[rr,1]*prob_dead[rr,cc]
    }
  }
# Never Seen
tempo3<-marray[[3]]
sumdo<-NULL
for(rr in 1:c(nocc-1)){
  sumdo<-0
  for(cc in c(rr+1):nocc) sumdo<-sumdo+prob_dead[rr,cc]+prob_alive[rr,cc]
  exp_ns[rr]<-tempo1[rr,1]*(1-sumdo)
 } 
expected<-list(expected_live=exp_alive,expected_dead=exp_dead,expected_neverseen=exp_ns)
####################### Calculate Cell Likelihoods ###############################
#alive
tempo1<-marray[[1]]
tempo2<-marray[[2]]
temp_ll1<-matrix(0,nrow=nrow(tempo1),ncol=nocc)
temp_ll2<-matrix(0,nrow=nrow(tempo2),ncol=nocc)
for(rr in 1:c(nocc-1)){
  first<-liveschedule[rr,2]
  last<-liveschedule[rr,3]
  for(cc in first:last){
    if(tempo1[rr,cc]>0) temp_ll1[rr,cc]<-tempo1[rr,cc]*log(max(constant,prob_alive[rr,cc]))
  }
}
for(rr in 1:c(nocc-1)){
  first<-deadschedule[rr,2]
  last<-deadschedule[rr,3]
  for(cc in first:last){
      if(tempo2[rr,cc]>0) temp_ll2[rr,cc]<-tempo2[rr,cc]*log(max(constant,prob_dead[rr,cc]))
  }
}
# Never Seen
tempo3<-marray[[3]]
sumdo<-NULL
temp_ll3<-rep(0,nocc-1)
for(rr in 1:c(nocc-1)){
  sumdo<-0
  for(cc in c(rr+1):nocc) sumdo<-sumdo+prob_dead[rr,cc]+prob_alive[rr,cc]
  if(tempo3[rr]>0) temp_ll3[rr]<-tempo3[rr]*log(1-sumdo)
} 
cell_likelihoods<-list(select_cell_live=temp_ll1,select_cell_dead=temp_ll2,cell_never_seen=temp_ll3)
######## Diagnostics #################################
deadcnt<-0
livecnt<-0
#alive
tempo1<-marray[[1]]
tempo2<-marray[[2]]
tempo3<-marray[[3]]
chi_alive<-matrix(0,nrow=nrow(tempo1),ncol=nocc)
pear_alive<-matrix(0,nrow=nrow(tempo1),ncol=nocc)
chi_dead<-matrix(0,nrow=nrow(tempo2),ncol=nocc)
pear_dead<-matrix(0,nrow=nrow(tempo2),ncol=nocc)
chi_ns<-rep(0,nocc-1)
pear_ns<-rep(0,nocc-1)
deadcnt<-0
livecnt<-0
nscnt<-0
for(rr in 1:c(nocc-1)){
  first<-liveschedule[rr,2]
  last<-liveschedule[rr,3]
  for(cc in first:last){
    livecnt<-livecnt+1
       if(exp_alive[rr,cc]>0){
         chi_alive[rr,cc]<-(tempo1[rr,cc]-exp_alive[rr,cc])^2/exp_alive[rr,cc]
         pear_alive[rr,cc]<-(tempo1[rr,cc]-exp_alive[rr,cc])/exp_alive[rr,cc]
       }
  }
}
for(rr in 1:c(nocc-1)){
  first<-deadschedule[rr,2]
  last<-deadschedule[rr,3]
  for(cc in first:last){
    deadcnt<-deadcnt+1
      if(exp_dead[rr,cc]>0){
        chi_dead[rr,cc]<-(tempo2[rr,cc]-exp_dead[rr,cc])^2/exp_dead[rr,cc]
        pear_dead[rr,cc]<-(tempo2[rr,cc]-exp_dead[rr,cc])/exp_dead[rr,cc]
      }
  }
}
# Never Seen
tempo3<-marray[[3]]
sumdo<-NULL
temp_ll3<-rep(0,nocc-1)
for(rr in 1:c(nocc-1)){
  sumdo<-0
  for(cc in c(rr+1):nocc) sumdo<-sumdo+prob_dead[rr,cc]+prob_alive[rr,cc]
  nscnt<-nscnt+1
     if(exp_ns[rr]>0){
      chi_ns[rr]<-(tempo3[rr]-exp_ns[rr])^2/exp_ns[rr]
      pear_ns[rr]<-(tempo3[rr]-exp_ns[rr])/exp_ns[rr]
    }
   
}
chisquare<-list(chi_live=chi_alive,chi_dead=chi_dead,chi_ns=chi_ns)
pearson<-list(pear_live=pear_alive,pear_dead=pear_dead,pear_ns=pear_ns)    
#### Test Statistics #################################
chi<-sum(chi_alive)+sum(chi_dead)+sum(chi_ns)
K<-length(results$par)
unpdf<-(deadcnt+livecnt+nscnt)-(nocc-1)-K
AIC<-2*results$value+2*K
sign<-1-pchisq(chi,unpdf)
N<-sum(marray[[1]][,1])
AICc<-AIC+(2*K*(K+1))/(sum(N)-K-1)
if(results$convergence==0) mess<-"Successful."
if(results$convergence==1) mess<-"Maximum Iterations Reached."
if(results$convergence==51) mess<-results$message
if(results$convergence==53) mess<-results$message
ans<-NULL
ans$statistics<-matrix(NA,7L,1L)
ans$statistics<-rbind(round(results$value*-1,3),round(K,0),round(AIC,2),round(AICc,2)
                      ,round(sum(N,0)),round(chi,2),round(unpdf,0))   
dimnames(ans$statistics)<-list(c("Log-Likelihood","K","AIC","AICc","Eff. Sample Size","Chi-square",
                                 "df"))
dimnames(ans$statistics)[[2]][1]<-list(c("Value"))
ans$model_convergence<-mess
ans$parameter_estimates<-parmtable
temps<-as.data.frame(Design)
names(temps)<-c("Occasion","Fdesign","Mdesign","Pdesign")
ans$design<-temps
ans$m_array<-marray
ans$expected<-expected
ans$cell_likelihoods<-cell_likelihoods
ans$chisquare<-chisquare
ans$pearson<-pearson
ans$prob_live<-prob_alive
ans$prob_dead<-prob_dead
ans$whichlivecells<-liveschedule
ans$whichdeadcells<-deadschedule
ans$type="hightower"
return(ans)
}
if(class(results)=="try-error") stop("Fit Failed.")
}