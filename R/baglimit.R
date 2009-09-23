baglimit<-function(intdir=NULL, estdir=NULL,species=NULL, state=NULL, mode=NULL, wave=NULL, styr=NULL, endyr=NULL,bag=NULL){
    if(is.null(intdir)) stop("Need main directory location of intercept files.")
    if(is.null(estdir)) stop("Need main directory location of catch and effort files.")
    if(is.null(species)) stop("Need NODC code for species.")
    if(is.null(state)) {warning("No state code was specified. Data from all states will be used."); state<-c(1:56,72,78)}
    if(is.null(mode)) {warning("No mode code was specified. Data for all modes will be used.");mode<-c(1,2,3,4,5,6,7,8,9)}
    if(is.null(wave)) {warning("No wave code was specified. Data for all wave will be used.");wave<-c(1,2,3,4,5,6)}
    if(is.null(bag)) stop("Bag limit value was not specified.")
    if(is.null(styr)) stop("Starting year is missing.")
    if(is.null(endyr)) stop("Ending year is missing.")

    if(length(grep("/",intdir))==1){
        din<-ifelse(substr(intdir,nchar(intdir),nchar(intdir)) %in% c("/"),
          c(paste(intdir,"int",sep="")),c(paste(intdir,"/int",sep="")))
     }
    if(length(grep("\\\\",intdir))==1){
        din<-ifelse(substr(intdir,nchar(intdir),nchar(intdir)) %in% c("\\"),
          c(paste(intdir,"int",sep="")),c(paste(intdir,"\\int",sep="")))
     }

     if(length(grep("/",estdir))==1){
        dest<-ifelse(substr(estdir,nchar(estdir),nchar(estdir)) %in% c("/"),
          c(paste(estdir,"est",sep="")),c(paste(estdir,"/est",sep="")))
     }
    if(length(grep("\\\\",estdir))==1){
        dest<-ifelse(substr(estdir,nchar(estdir),nchar(estdir)) %in% c("\\"),
          c(paste(estdir,"est",sep="")),c(paste(estdir,"\\est",sep="")))
     }
 if(!all(mode %in% c(3,4,5,6,7))) stop ("Mode not valid.")
 if(!all(wave %in% c(1:6))) stop ("Wave not valid.")
if(any(state==121) & any(state==122)) stop("Use state 12 to combine Florida coasts.")
if(any(state==12) & any(state %in% c(121,122))) stop("You are mixing codes for Florida.")

if(styr>=1982 & endyr<=2004){
   if(!any(state %in% c(1,12,13,22,28,37,45,121,122)) & any(mode %in% c(4,5)))
       warning("Mode codes 4 or 5 not available for a state during the selected years.")
  }
if(styr>=2005 | endyr>=2005){
   if(any(mode %in% c(6)))
     warning ("Mode code 6 for combined party/charter boats not available after 2005.")
  }
if(!all(state %in% c(1,2,4,5,6,8,9,10,11,12,13,15:42,44:51,53:56,90,121,122))) stop ("Invalid state code.")
ID_CODE<-NULL;FSHINSP<-NULL;ST<-NULL;LEADER<-NULL;NUM_FISH<-NULL;
MODE_FX<-NULL;CNTRBTRS<-NULL;CATCH<-NULL;NUMRTRIP<-NULL
species<-as.character(species)
	flag<-0;outresults<-NULL
	for(yr in styr:endyr){
 	   for (j in 1:as.numeric(length(wave))){ 
      	  wv<-wave[j]   
              t3<-read.csv(paste(din,yr,"/","I3_",yr,wv,".csv",sep=""))
              t3$ID_CODE<-as.character(t3$ID_CODE) 
              t3$SP_CODE<-as.character(t3$SP_CODE)  

       if(yr==1992 & wv==1){
          t3$WAVE<-ifelse(t3$SUB_REG>=6 & t3$MONTH==1,1.1,
	          ifelse(t3$SUB_REG>=6 & t3$MONTH==2,1.2,t3$WAVE))
          t3$del<-ifelse(t3$SUB_REG>=6 & t3$MONTH>2,1,0)
          t3<-t3[t3$del==0,]
        }
       if(yr==1988){
          if(wv==4){
             t3$WAVE<-ifelse(t3$SUB_REG>=6 &t3$MONTH==7,4.1,
                    ifelse(t3$SUB_REG>=6 & t3$MONTH==8,4.2,t3$WAVE))
             t3$del<-ifelse(t3$SUB_REG>=6 & (t3$MONTH<7|t3$MONTH>8),1,0)
             t3<-t3[t3$del==0,]
           }
          if(wv==5){
             t3$WAVE<-ifelse(t3$SUB_REG>=6 &t3$MONTH==9,5.1,
                  ifelse(t3$SUB_REG>=6 & t3$MONTH==10,5.2,t3$WAVE))
             t3$del<-ifelse(t3$SUB_REG>=6 & (t3$MONTH<9|t3$MONTH>10),1,0)
             t3<-t3[t3$del==0,]
           }
          if(wv==6){
             t3$WAVE<-ifelse(t3$SUB_REG>=6 & t3$MONTH==11,6.1,
                ifelse(t3$SUB_REG>=6 & t3$MONTH==12,6.2,t3$WAVE))
             t3$del<-ifelse(t3$SUB_REG>=6 & (t3$MONTH<11|t3$MONTH>12),1,0)
             t3<-t3[t3$del==0,]
           }
       }
              t3$MODE_FX<-ifelse(t3$YEAR==2004 & t3$SUB_REG<6 & 
                    t3$MODE_F==6,6,t3$MODE_FX)
              t3$MODE_FX<-ifelse(t3$YEAR==2004 & t3$SUB_REG<6 & 
                    t3$MODE_F==7,6,t3$MODE_FX)
              t3$MODE_FX<-ifelse(t3$YEAR==2005 & t3$SUB_REG<6 & 
                    t3$MODE_F==6,4,t3$MODE_FX)
               t3$MODE_FX<-ifelse(t3$YEAR==2005 & t3$SUB_REG<6 & 
                    t3$MODE_F==7,5,t3$MODE_FX)
              t3$MODE_FX<-ifelse(t3$YEAR==2006 & t3$SUB_REG<6 & 
                    t3$MODE_F==6,4,t3$MODE_FX)
               t3$MODE_FX<-ifelse(t3$YEAR==2006 & t3$SUB_REG<6 & 
                    t3$MODE_F==7,5,t3$MODE_FX)
               t3$MODE_FX<-ifelse(t3$MODE_FX<3,3,t3$MODE_FX)

        	  if(any(state %in% c(121,122))){
                      d3<-t3[t3$ST!=12,]
                      if(any(state==121)) d4<-t3[t3$ST==12 & t3$SUB_REG!=7,]
                      if(any(state==122)) d4<-t3[t3$ST==12 & t3$SUB_REG!=6,]
                      t3<-as.data.frame(rbind(d3,d4))
                      state1<-state[state!=121|state!=122]
                      state1<-c(state1,12)
                  }
            if(!any(state %in% c(121,122))) state1<-state
        	  t3<-t3[t3$SP_CODE==species & t3$ST %in% c(state1) & 
              t3$MODE_FX %in% c(mode) & t3$AREA_X %in% c(1:5),]
              t3$cnt<-ifelse(is.na(t3$LNGTH),NA,1)
              if(length(t3$ID_CODE>0)){	
              t3<-aggregate(t3$cnt,list(t3$ID_CODE,t3$FSHINSP),sum,na.rm=T)
              t3[,2]<-as.numeric(as.character(t3[,2]))
              t3[,2]<-ifelse(t3[,3]>t3[,2],t3[,3],t3[,2])
              t3<-t3[,1:2]
              names(t3)<-c("ID_CODE","FSHINSP")
              }
              t3<-t3[duplicated(t3$ID_CODE)==FALSE,]
              t3<-subset(t3,select=c(ID_CODE,FSHINSP))
 
             if(length(t3$ID_CODE)==0|is.na(t3$ID_CODE[1])){
               warning(paste("Species not found in wave ",wv))  
               next
            }

        	  t2<-read.csv(paste(din,yr,"/","I2_",yr,wv,".csv",sep=""))
	  	  t2$ID_CODE<-as.character(t2$ID_CODE) 
        	  t2$SP_CODE<-as.character(t2$SP_CODE)
 
            if(length(t2$ID_CODE)==0|is.na(t2$ID_CODE[1])){
               warning(paste("Species not found in wave ",wv))  
               next
            }


     t4<-read.csv(paste(din,yr,"/","I4_",yr,wv,".csv",sep=""))
           t4$ID_CODE<-as.character(t4$ID_CODE) 
           t4$LEADER<-as.character(t4$LEADER) 
           t4<-subset(t4,ST>0,select=c(ID_CODE,LEADER))
           t4<-t4[order(t4$LEADER),]

        #Adjust type2 for leader code
	     adj_t2<-merge(t2,t4,by.x="ID_CODE",by.y="ID_CODE",all.x=T)
           adj_t2$oldid<-ifelse(!is.na(adj_t2$ID_CODE) & !is.na(adj_t2$LEADER),
                  adj_t2$ID_CODE,NA)
           adj_t2$ID_CODE<-ifelse(!is.na(adj_t2$ID_CODE) & !is.na(adj_t2$LEADER),
               adj_t2$LEADER,adj_t2$ID_CODE)
           rm(t2,t4)
           t2<-adj_t2[order(adj_t2$ID_CODE),]

             if(yr==1992 & wv==1){
                 t2$WAVE<-ifelse(t2$SUB_REG>=6 & t2$MONTH==1,1.1,
	                ifelse(t2$SUB_REG>=6 & t2$MONTH==2,1.2,t2$WAVE))
                 t2$del<-ifelse(t2$SUB_REG>=6 & t2$MONTH>2,1,0)
                 t2<-t2[t2$del==0,]
                }
             if(yr==1988){
                 if(wv==4){
                    t2$WAVE<-ifelse(t2$SUB_REG>=6 &t2$MONTH==7,4.1,
                      ifelse(t2$SUB_REG>=6 & t2$MONTH==8,4.2,t2$WAVE))
                      t2$del<-ifelse(t2$SUB_REG>=6 & (t2$MONTH<7|t2$MONTH>8),1,0)
                      t2<-t2[t2$del==0,]
                  }
                if(wv==5){
                    t2$WAVE<-ifelse(t2$SUB_REG>=6 &t2$MONTH==9,5.1,
                      ifelse(t2$SUB_REG>=6 & t2$MONTH==10,5.2,t2$WAVE))
                    t2$del<-ifelse(t2$SUB_REG>=6 & (t2$MONTH<9|t2$MONTH>10),1,0)
                    t2<-t2[t2$del==0,]
                 }
               if(wv==6){
                   t2$WAVE<-ifelse(t2$SUB_REG>=6 &t2$MONTH==11,6.1,
                         ifelse(t2$SUB_REG>=6 & t2$MONTH==12,6.2,t2$WAVE))
                   t2$del<-ifelse(t2$SUB_REG>=6 & (t2$MONTH<11|t2$MONTH>12),1,0)
                  t2<-t2[t2$del==0,]
                }
             }
               t2$MODE_FX<-ifelse(t2$YEAR==2004 & t2$SUB_REG<6 & 
              t2$MODE_F==6,6,t2$MODE_FX)
              t2$MODE_FX<-ifelse(t2$YEAR==2004 & t2$SUB_REG<6 & 
                    t2$MODE_F==7,6,t2$MODE_FX)
              t2$MODE_FX<-ifelse(t2$YEAR==2005 & t2$SUB_REG<6 & 
                    t2$MODE_F==6,4,t2$MODE_FX)
               t2$MODE_FX<-ifelse(t2$YEAR==2005 & t2$SUB_REG<6 & 
                    t2$MODE_F==7,5,t2$MODE_FX)
              t2$MODE_FX<-ifelse(t2$YEAR==2006 & t2$SUB_REG<6 & 
                    t2$MODE_F==6,4,t2$MODE_FX)
               t2$MODE_FX<-ifelse(t2$YEAR==2006 & t2$SUB_REG<6 & 
                    t2$MODE_F==7,5,t2$MODE_FX)
               t2$MODE_FX<-ifelse(t2$MODE_FX<3,3,t2$MODE_FX)

                 if(any(state %in% c(121,122))){
                      d3<-t2[t2$ST!=12,]
                      if(any(state==121)) d4<-t2[t2$ST==12 & t2$SUB_REG!=7,]
                      if(any(state==122)) d4<-t2[t2$ST==12 & t2$SUB_REG!=6,]
                      t2<-as.data.frame(rbind(d3,d4))
                      state1<-state[state!=121|state!=122]
                      state1<-c(state1,12)
                  }
            if(!any(state %in% c(121,122))) state1<-state
        	
         t2<-t2[t2$SP_CODE==species & t2$ST %in% c(state1) & t2$MODE_FX %in% c(mode) & 
                  t2$DISPO %in% c(3:7) &  t2$AREA_X %in% c(1:5),]
         t2<-t2[t2$CNTRBTRS!=0,]
          t2<-t2[duplicated(t2$ID_CODE,t2$DISPO)==FALSE,]
           if(length(t2$ID_CODE)>0){
              t2<-aggregate(t2$NUM_FISH,list(t2$ID_CODE),sum,na.rm=T)
               names(t2)<-c("ID_CODE","NUM_FISH")
                 t2$ID_CODE<-as.character(t2$ID_CODE)
             }
          t2<-subset(t2,select=c(ID_CODE,NUM_FISH))
          m1<-merge(t2,t3,by.x=c("ID_CODE"),by.y=c("ID_CODE"),sort=T,all.x=T,all.y=T)  
          m1<-m1[!is.na(m1$ID_CODE),] 
 
          if(length(m1$ID_CODE)==0){
               warning(paste("Species not found in wave ",wv))  
               next
            }
       t1<-read.csv(paste(din,yr,"/","I1_",yr,wv,".csv",sep=""))
            t1$ID_CODE<-as.character(t1$ID_CODE)
       if(yr==1992 & wv==1){
          t1$WAVE<-ifelse(t1$SUB_REG>=6 & t1$MONTH==1,1.1,
	          ifelse(t1$SUB_REG>=6 & t1$MONTH==2,1.2,t1$WAVE))
          t1$del<-ifelse(t1$SUB_REG>=6 & t1$MONTH>2,1,0)
          t1<-t1[t1$del==0,]
        }
       if(yr==1988){
          if(wv==4){
             t1$WAVE<-ifelse(t1$SUB_REG>=6 &t1$MONTH==7,4.1,
                    ifelse(t1$SUB_REG>=6 & t1$MONTH==8,4.2,t1$WAVE))
             t1$del<-ifelse(t1$SUB_REG>=6 & (t1$MONTH<7|t1$MONTH>8),1,0)
             t1<-t1[t1$del==0,]
           }
          if(wv==5){
             t1$WAVE<-ifelse(t1$SUB_REG>=6 &t1$MONTH==9,5.1,
                  ifelse(t1$SUB_REG>=6 & t1$MONTH==10,5.2,t1$WAVE))
             t1$del<-ifelse(t1$SUB_REG>=6 & (t1$MONTH<9|t1$MONTH>10),1,0)
             t1<-t1[t1$del==0,]
           }
          if(wv==6){
             t1$WAVE<-ifelse(t1$SUB_REG>=6 &t1$MONTH==11,6.1,
                ifelse(t1$SUB_REG>=6 & t1$MONTH==12,6.2,t1$WAVE))
             t1$del<-ifelse(t1$SUB_REG>=6 & (t1$MONTH<11|t1$MONTH>12),1,0)
             t1<-t1[t1$del==0,]
           }
       }
               t1$MODE_FX<-ifelse(t1$YEAR==2004 & t1$SUB_REG<6 & 
                    t1$MODE_F==6,6,t1$MODE_FX)
               t1$MODE_FX<-ifelse(t1$YEAR==2004 & t1$SUB_REG<6 & 
                    t1$MODE_F==7,6,t1$MODE_FX)
              t1$MODE_FX<-ifelse(t1$YEAR==2005 & t1$SUB_REG<6 & 
                    t1$MODE_F==6,4,t1$MODE_FX)
               t1$MODE_FX<-ifelse(t1$YEAR==2005 & t1$SUB_REG<6 & 
                    t1$MODE_F==7,5,t1$MODE_FX)
              t1$MODE_FX<-ifelse(t1$YEAR==2006 & t1$SUB_REG<6 & 
                    t1$MODE_F==6,4,t1$MODE_FX)
               t1$MODE_FX<-ifelse(t1$YEAR==2006 & t1$SUB_REG<6 & 
                    t1$MODE_F==7,5,t1$MODE_FX)
               t1$CNTRBTRS<-ifelse(t1$CNTRBTRS %in% c(0,88,99,NA) & t1$NUM_TYP4==0,1,
                       t1$CNTRBTRS)
               t1$CNTRBTRS<-ifelse(t1$CNTRBTRS %in% c(0,88,99,NA) & t1$NUM_TYP4!=0,0,
                         t1$CNTRBTRS) 
                t1<-t1[t1$CNTRBTRS!=0,]
                if(any(state %in% c(121,122))){
                      d3<-t1[t1$ST!=12,]
                      if(any(state==121)) d4<-t1[t1$ST==12 & t1$SUB_REG!=7,]
                      if(any(state==122)) d4<-t1[t1$ST==12 & t1$SUB_REG!=6,]
                      t1<-as.data.frame(rbind(d3,d4))
                      state1<-state[state!=121|state!=122]
                      state1<-c(state1,12)
                  }
            if(!any(state %in% c(121,122))) state1<-state
       	  t1<-t1[t1$ST %in% c(state1) & t1$MODE_FX %in% c(mode) &  t1$AREA_X %in% c(1:5),]
        	  t1<-t1[order(t1$ID_CODE),names(t1)!="NUM_FISH"]
              t1<-t1[duplicated(t1$ID_CODE)==FALSE,]
                      
        if(length(m1$ID_CODE)>0){
          alldata<-merge(m1,t1,by.x="ID_CODE",by.y="ID_CODE",all.x=T,all.y=T)
          alldata<-alldata[duplicated(t1$ID_CODE)==FALSE,]
            for(k in 1:ncol(alldata)){
               if(is.numeric(alldata[,k])) alldata[is.na(alldata[,k]),k] <- 0
              } 
             alldata$CATCH<-alldata$NUM_FISH+alldata$FSHINSP
             alldata<-alldata[alldata$CNTRBTRS!=0,]
            #Make a record for each person in group with cntrbrs >1
             getdata<-alldata[alldata$CNTRBTRS>1,]
           if(length(getdata$ID_CODE)>0){
             getdata$CATCH<-ifelse(is.nan(getdata$CATCH/getdata$CNTRBTRS),0,getdata$CATCH/getdata$CNTRBTRS)
             getdata<-subset(getdata,select=c(ID_CODE,ST,MODE_FX,CNTRBTRS,CATCH))
              newdata<-NULL
              for(i in 1:length(getdata$ID_CODE)){
                 loop<-getdata[i,4]
                 for(j in 1:loop){
                   newdata<-rbind(newdata,getdata[i,])
                  }
               }  
              newdata<-newdata[,c(1,2,3,4,5)]
              alldata<-subset(alldata,select=c(ID_CODE,ST,MODE_FX,CNTRBTRS,CATCH)) 
              alldata<-alldata[alldata$CNTRBTRS==1,]
              alldata<-rbind(alldata,newdata)
            }  
            alldata<-subset(alldata,select=c(ID_CODE,ST,MODE_FX,CNTRBTRS,CATCH))     
        	ti<-as.data.frame(table(alldata$ST,alldata$MODE_FX))
        	  names(ti)<-c("ST","MODE_FX","IF")     
  
       	  cf<-as.data.frame(table(alldata$ST,alldata$MODE_FX,alldata$CATCH))
        	  names(cf)<-c("ST","MODE_FX","CATCH","CF")
              cf<-cf[order(cf$MODE_FX,cf$CATCH),]
              cf<-cf[cf$CATCH!="0",]
        	  fd<-merge(cf,ti,by.x=c("ST","MODE_FX"),by.y=c("ST","MODE_FX"))
        	  fd$REL<-ifelse(is.nan(fd$CF/fd$IF),0,fd$CF/fd$IF)

         	ests<-read.csv(paste(dest,yr,"/","AG_",yr,wv,".csv",sep="")) 
         	ests$SPCODE<-as.character(ests$SP_CODE)
            ests<-ests[ests$POOL_FLG==1 & ests$OUTFLG==1 & ests$EX_FLG==1,] 
            ests$WAVE<-ifelse(ests$WAVE %in% c(4.1,4.2),4,ifelse(
                              ests$WAVE %in% c(5.1,5.2),5,ifelse(
                              ests$WAVE %in% c(6.1,6.2),6,ests$WAVE)))                                                                                                                                                                           
            ests$MODE_FX<-ifelse(ests$MODE_FX<3,3,ests$MODE_FX)
               
          if(any(state %in% c(121,122))){
                      d3<-ests[ests$ST!=12,]
                      if(any(state==121)) d4<-ests[ests$ST==12 & ests$SUB_REG!=7,]
                      if(any(state==122)) d4<-ests[ests$ST==12 & ests$SUB_REG!=6,]
                      ests<-as.data.frame(rbind(d3,d4))
                      state1<-state[state!=121|state!=122]
                      state1<-c(state1,12)
                  }
            if(!any(state %in% c(121,122))) state1<-state
	   	ests<-ests[ests$ST %in% c(state1) & ests$MODE_FX %in% c(mode)
                  & ests$AREA_X %in% c(1:5),]

         	ests<-unique(subset(ests,select=c(ST, MODE_FX, NUMRTRIP)))
 
   	      ffd<-merge(fd,ests,by.x=c("ST","MODE_FX"),by.y=c("ST","MODE_FX"))
        	  ffd$TRIPS<-ifelse(is.nan(ffd$REL*ffd$NUMRTRIP),0,round(ffd$REL*ffd$NUMRTRIP,1))
            if(length(ffd$ST)>0){
       	f1<-aggregate(ffd$TRIPS,list(ffd$CATCH),sum);names(f1)<-c("CATCH","TRIPS") 
        	f2<-aggregate(ffd$CF,list(ffd$CATCH),sum);names(f2)<-c("CATCH","CF") 
        	final<-merge(f1,f2,by.x="CATCH",by.y="CATCH")
          		 final$YEAR<-yr;final$WAVE<-wv
                  flag<-flag+1 
              }
        	if(flag==1) outresults<-final
       	if(flag>1) outresults<-rbind(outresults,final)
       }
      }
    }
if(is.null(outresults)) stop ("No data found.")

       outres<-aggregate(outresults[,2:3],list(outresults$CATCH),sum)     
       names(outres)<-c("CATCH","TRIPS","INTERCEPTS")
             outres$CATCH<-as.numeric(levels(outres$CATCH))[as.integer(outres$CATCH)] 
              ccat<-as.data.frame(1:max(ceiling(outres$CATCH)));names(ccat)<-c("CATCH")
              outres<-merge(ccat,outres,by.x="CATCH",by.y="CATCH",all.x=T,all.y=T)
              outres$TRIPS<-ifelse(is.na(outres$TRIPS),0,outres$TRIPS)
              outres$INTERCEPTS<-ifelse(is.na(outres$INTERCEPTS),0,outres$INTERCEPTS)
               outres$PERCENT<-round(outres$INTERCEPTS/sum(outres$INTERCEPTS)*100,3) 
              outres$CUMPER<-round(cumsum(outres$PERCENT),1)
       outres$ORIG.HARVEST<-ifelse(is.nan(outres$CATCH*outres$TRIPS),0,round(outres$CATCH*outres$TRIPS,1))
  
       pred<-NULL
       outpt<-list(outres);names(outpt)<-c("No Bag Limit")
       for(i in 1:as.numeric(length(bag))){
           bag<-sort(bag)
           if(bag[i]<=max(outres$CATCH)){
              sums<-apply(outres[outres[,1]>=bag[i],c(2,3,4)],2,sum)
              bg<-outres[,1:5];bg[which(bg$CATCH==bag[i]),2]<-sums[1];bg[which(bg$CATCH==bag[i]),3]<-sums[2]
              bg[which(bg$CATCH==bag[i]),4]<-sums[3]
              bg[as.numeric(which(bg$CATCH==bag[i])+1):as.numeric(length(outres$CATCH)),2:4]<-0
              bg$CUMPER<-round(cumsum(bg$PERCENT),1)
              bg$HARVEST<-ifelse(is.nan(bg$CATCH*bg$TRIPS),0,round(bg$CATCH*bg$TRIPS,1))
              pred[i]<-round((sum(bg$HARVEST)-sum(outres$ORIG.HARVEST))/sum(outres$ORIG.HARVEST)*100,2)
              outpt[[i+1]]<-bg;names(outpt)[i+1]<-c(paste("Bag Limit:",bag[i],sep=""))
           }
          if(bag[i]>max(outres$CATCH)){
           pred[i]<-0
           out1<-outres;names(out1)[6]<-c("HARVEST")
           outpt[[i+1]]<-out1
           names(outpt)[i+1]<-c(paste("Bag Limit:",bag[i],sep=""))
         }
      }
      one<-as.data.frame(cbind(bag,abs(pred)));names(one)<-c("Bag Limit","% Reduction")
      outpt$Results<-one
              details<-NULL
             details[1]<-c(paste("Species: ",species,sep=""))
             details[2]<-c(paste("State: ",paste(state,collapse=","),sep=""))
             details[3]<-c(paste("Start and End Years: ",styr,"-",endyr,sep=""))
             details[4]<-c(paste("Wave: ",paste(wave,collapse=","),sep=""))
             details[5]<-c(paste("Mode: ",paste(mode,collapse=","),sep=""))
             outpt$Details<-details
       return(outpt)
}

#baglimit(intdir="C:/Temp", estdir="C:/Temp",species=8835250101, state=36, mode=c(7), wave=c(8), styr=2007, endyr=2007,bag=c(5,10))

