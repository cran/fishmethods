
fm_checkdesign<-function(occasions=NULL,design=NULL,type="F"){
#################### Setup Fdesign ####################################
nocc<-occasions
Design<-matrix(0,nrow=c(nocc-1),2)
pardef<-strsplit(design,"*",fixed=TRUE)
yes1<-0;yes2<-0;seqs<-NULL;subcells<-NULL
if(type %in% c("F","M")){
  Design[,1]<-1:c(nocc-1)
  pat<-as.character(nocc)
  for(t in 1:length(pardef)){
      tempo<-pardef[[t]]
      if(length(tempo)==1) {
        eval(parse(text=paste("seqs<-",tempo[1],sep="")))
        if(any(seqs==pat)) yes1<-1 
        if(yes1==1) stop(paste("The occasions for F or M range from 1 to",nocc-1,sep=" "))
      }
      if(length(tempo)>1){
        eval(parse(text=paste("seqs<-",tempo[1],sep="")))
        eval(parse(text=paste("subcells<-",tempo[2],sep="")))
        if(any(seqs==pat)) yes1<-1
        if(any(subcells==pat)) yes2<-1
        if(yes1==1 || yes2==1) stop(paste("The occasions for F or M range from 1 to",nocc-1,sep=" "))
      }
  }
  pat<-1
  yes1<-0
  for(t in 1:length(pardef)){
    tempo<-pardef[[t]]
    if(length(tempo)==1) {
      eval(parse(text=paste("seqs<-",tempo[1],sep="")))
      if(any(seqs==1)) yes1<-1 
    }
    if(length(tempo)>1){
      eval(parse(text=paste("seqs<-",tempo[1],sep="")))
      eval(parse(text=paste("subcells<-",tempo[2],sep="")))
      if(any(seqs==pat)) yes1<-1
      if(any(subcells==pat)) yes1<-1
    }
  }
  if(yes1==0) stop(paste("Missing occasion 1 for F or M"))
}  
if(type %in% c("P")){
  yes1<-0
  Design[,1]<-2:c(nocc)
  pat<-1
  for(t in 1:length(pardef)){
    tempo<-pardef[[t]]
    if(length(tempo)==1){
      eval(parse(text=paste("seqs<-",tempo[1],sep="")))
      if(any(seqs==1)) yes1<-1 
      if(yes1==1) stop(paste("The occasions for P range from 2 to",nocc,sep=" "))
    }  
    if(length(tempo)>1){
      eval(parse(text=paste("seqs<-",tempo[1],sep="")))
      eval(parse(text=paste("subcells<-",tempo[2],sep="")))
      if(any(seqs==pat)) yes1<-1
      if(any(subcells==pat)) yes1<-1
      if(yes1==1) stop(paste("The occasions for P range from 2 to",nocc,sep=" "))
    }
  }
}  
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
if(Design[1,2]>0){
  for(t in 1:length(Design[,1])){
  if(Design[t,2]==0) Design[t,2]<-Design[t-1,2]
  }
}
Design<-as.data.frame(Design)
names(Design)<-c("Occasions","Parameter")
return(Design)
}

