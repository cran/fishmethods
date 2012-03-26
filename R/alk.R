
alk<-function(age=NULL,size=NULL,binsize=NULL,type=1){
         if(is.null(age)) stop ("age vector does not exist")
         if(!is.numeric(age)) stop ("age vector is not numeric")
         if(is.null(size)) stop ("size vector does not exist")
         if(!is.numeric(size)) stop ("size vector is not numeric")
         d<-NULL;outs<-NULL;ll<-NULL;ul<-NULL;la<-NULL;ua<-NULL;lenlist<-NULL
         agelist<-NULL
         d<-as.data.frame(cbind(age,size))
         d<-d[!is.na(d$age),]
	    d<-d[!is.na(d$size),]
        
         d$lenclass<-trunc(size/binsize)*binsize+binsize/2
         ll<-min(d$lenclass)
         ul<-max(d$lenclass)
         la<-min(d$age)
         ua<-max(d$age)
         lenlist<-seq(ll,ul,by=binsize)
         agelist<-seq(la,ua,by=1)
	    agelen<-matrix(0,nrow=length(lenlist),ncol=length(agelist))
          for(loops in 1:as.numeric(length(d$age))){
              for(a in 1:length(agelist)){
                for(b in 1:length(lenlist)){
                   if(d$age[loops]==agelist[a] & d$lenclass[loops]==lenlist[b]){
                       agelen[b,a]<-agelen[b,a]+1
                   }
                }
             }
           }
 	if(type==1){	
       outs<-cbind(lenlist,rowSums(agelen),as.data.frame(agelen))
       names(outs)<-c("len","nl",paste("A",agelist,sep=""))
     }
     if(type==2){
       outs<-cbind(lenlist,rowSums(agelen),as.data.frame(agelen)/rowSums(agelen))
       names(outs)<-c("len","nl",paste("A",agelist,sep=""))
       for(i in 3:ncol(outs)){
         for(j in 1:nrow(outs)){
             if(is.nan(outs[j,i])) outs[j,i]<-0
		}
	  }
     }  
     return(outs)
 }

