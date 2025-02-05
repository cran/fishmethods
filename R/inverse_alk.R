inverse_alk<-function(alk1=NULL,lf1=NULL,lf2=NULL,toler=0.000001,max.iter=10000){
  if(is.null(alk1)) stop("alk1 is missing")
  if(is.null(lf2)) stop("lf2 is missing")
  if(!is.matrix(alk1)) stop("alk1 must be a matrix")
  if(dim(alk1)[2]<dim(alk1)[1]) stop("Number of size columns must be >= the number of age rows")
  if(!is.null(lf1)){
    if(length(lf1)!=ncol(alk1)) stop("Number of size intervals differs between alk1 matrix and lf1 vector")
      denom<-colSums(alk1)
    alk1<-round(t(lf1*(t(alk1)/denom)),0)
  }
  nages<-nrow(alk1)
  nlens<-ncol(alk1)
  n_plus_j_1<-colSums(alk1)
  n_i_plus_1<-rowSums(alk1)
  n_plus_plus_1<-sum(alk1)
  if(length(n_plus_j_1)!=nlens) stop("Number of size intervals differs between alk1 matrix and lf2 vector")
  ai1<-vector()
  for(i in 1:nages) ai1[i]<-n_i_plus_1[i]/n_plus_plus_1
  n_plus_j_2<-as.vector(lf2)
  n_plus_plus_2<-sum(n_plus_j_2)
  alk2<-matrix(0,nrow=nages,ncol=nlens)
  Pj_i<-matrix(0,nrow=nages,ncol=nlens)
  ai2<-vector()
  LL<-matrix(0,nrow=nages,ncol=nlens)
  milf2<-matrix(0,nrow=nages,ncol=nlens)
  milf1<-matrix(0,nrow=nages,ncol=nlens)
  
########################## Solve ######################################
  Lprime<-0
  ok <- 0
  fuzz <- toler
  cnter<-0
#----------Initial Guesses--------------------#
  for(i in 1:nages){
    for(j in 1:nlens){
      alk2[i,j]<-n_plus_j_2[j]*alk1[i,j]/n_plus_j_1[j]
    }
  }
  n_i_plus_2<-rowSums(alk2)
#-------Iterations-------------------------#
while (ok == 0 & cnter<max.iter) {
#M step
  for(i in 1:nages){
    for(j in 1:nlens){
     Pj_i[i,j]<-(alk1[i,j]+alk2[i,j])/(n_i_plus_1[i]+n_i_plus_2[i])
    }
  }
  ai2<-n_i_plus_2/n_plus_plus_2
  
  for(i in 1:nages){
    for(j in 1:nlens){
     milf2[i,j]<-Pj_i[i,j]*n_i_plus_2[i]
    }
  }
  
#E step
  for(i in 1:nages){
    for(j in 1:nlens){
      alk2[i,j]<-(milf2[i,j]*n_plus_j_2[j])/sum(milf2[,j])
    }
  }
  n_i_plus_2<-rowSums(alk2)
  
  LL<-0
  for(i in 1:nages){
    for(j in 1:nlens){
      if(Pj_i[i,j]>0){
       LL<-LL+(alk1[i,j]+alk2[i,j])*log(Pj_i[i,j])
      }
    }
    if(ai2[i]>0){
    LL<-LL+n_i_plus_2[i]*log(ai2[i])
    }
  }
  cnter<-cnter+1 
  diff<-abs(Lprime-LL)
  if(diff<=fuzz) ok<-1
  if(diff>fuzz) Lprime<-LL
}#while

for(i in 1:nages){
  for(j in 1:nlens){
    milf1[i,j]<-Pj_i[i,j]*n_i_plus_1[i]
  }
}
resids<-(alk1-milf1)/sqrt(milf1)
resids[is.nan(resids)]<-0

i_prop<-round(data.frame(prop=rowSums(alk2)/sum(alk2)),6)


output<-list(alk1=alk1,lf1=lf1,lf2=lf2,est_alk2=alk2,residuals=resids,
             est_age_prop_lf2=i_prop)

return(output)
}#inverse_alk function

