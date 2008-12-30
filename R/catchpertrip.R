##############################################################################
#
#     Program: catch-per-trip analysis
#       CONDUCT CATCH PER TRIP ANALYSIS FOLLOWING MRFSS HANDBOOK
#
#       Gary Nelson
#       Massachusetts Division of Marine Fisheries
#       Gloucester, MA 01930
#       gary.nelson@state.ma.us
#
#############################################################################

catchpertrip<-function(intdir=NULL, estdir=NULL,species=NULL, state=NULL, mode=NULL, wave=NULL, styr=NULL, endyr=NULL){
    if(is.null(intdir)) stop("Need main directory location of intercept files.")
    if(is.null(estdir)) stop("Need main directory location of catch and effort files.")
    if(is.null(species)) stop("Need NODC code for species.")
    if(is.null(state)) {warning("No state code was specified. Data from all states will be used."); state<-c(1:56,72,78)}
    if(is.null(mode)) {warning("No mode code was specified. Data for all modes will be used.");mode<-c(1,2,3,4,5,6,7,8,9)}
    if(is.null(wave)) {warning("No wave code was specified. Data for all wave will be used.");wave<-c(1,2,3,4,5,6)}
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

  species<-as.character(species)
	flag<-0
	for(yr in styr:endyr){
 	   for (j in 1:as.numeric(length(wave))){ 
      	  wv<-wave[j] 
              t3<-read.csv(paste(din,yr,"/","I3_",yr,wv,".csv",sep=""))
              t3$ID_CODE<-as.character(t3$ID_CODE) 
              t3$SP_CODE<-as.character(t3$SP_CODE) 
if(yr==1992 & wv==1){
          t3$WAVE<-ifelse(t3$SUB_REG>=6 & t3$MONTH==1,1.1,
	          ifelse(t3$SUB_REG>=6 & t3$MONTH==2,1.2,t3$WAVE))
          t2$del<-ifelse(t3$SUB_REG>=6 & t3$MONTH>2,1,0)
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


  t2<-read.csv(paste(din,yr,"/","I2_",yr,wv,".csv",sep=""))
	  	  t2$ID_CODE<-as.character(t2$ID_CODE) 
        	  t2$SP_CODE<-as.character(t2$SP_CODE)

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
                  t2$DISPO %in% c(1:7) &  t2$AREA_X %in% c(1:5),]
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
          alldata<-alldata[alldata$CNTRBTRS==1,]

            for(k in 1:ncol(alldata)){
               if(is.numeric(alldata[,k])) alldata[is.na(alldata[,k]),k] <- 0
              } 
            alldata$CATCH<-alldata$SUM_FISH+alldata$FSHINSP
            flag<-flag+1 
        	ti<-as.data.frame(table(alldata$ST,alldata$MODE_FX))
        	  names(ti)<-c("ST","MODE_FX","IF")     
       	cf<-as.data.frame(table(alldata$ST,alldata$MODE_FX,alldata$CATCH))
        	  names(cf)<-c("ST","MODE_FX","CATCH","CF")
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

   		ests<-ests[ests$SP_CODE==species & ests$ST %in% c(state1) & ests$MODE_FX %in% c(mode),]
            ests<-unique(subset(ests,select=c(ST, MODE_FX, NUMRTRIP)))
        	ffd<-merge(fd,ests,by.x=c("ST","MODE_FX"),by.y=c("ST","MODE_FX"))
        	  ffd$TRIPS<-ifelse(is.nan(ffd$REL*ffd$NUMRTRIP),0,round(ffd$REL*ffd$NUMRTRIP,1))
        	f1<-aggregate(ffd$TRIPS,list(ffd$CATCH),sum);names(f1)<-c("CATCH","TRIPS") 
        	f2<-aggregate(ffd$CF,list(ffd$CATCH),sum);names(f2)<-c("CATCH","CF") 
        	final<-merge(f1,f2,by.x="CATCH",by.y="CATCH")
          		 final$YEAR<-yr;final$WAVE<-wv
        	if(flag==1) outresults<-final
       	if(flag>1) outresults<-rbind(outresults,final)
       }
      }
    }
if(is.null(outresults)) stop ("No data found.")

       outres<-aggregate(outresults[,2:3],list(outresults$CATCH),sum)     
       names(outres)<-c("CATCH","TRIPS","CF")
             outres$CATCH<-as.numeric(levels(outres$CATCH))[as.integer(outres$CATCH)] 
              ccat<-as.data.frame(1:max(outres$CATCH));names(ccat)<-c("CATCH")
              outres<-merge(ccat,outres,by.x="CATCH",by.y="CATCH",all.x=T,all.y=T)
              outres$TRIPS<-ifelse(is.na(outres$TRIPS),0,outres$TRIPS)
              outres$CF<-ifelse(is.na(outres$CF),0,outres$CF)
              outres$PROP_TRIP<-round(outres$TRIPS/sum(outres$TRIPS),3)
           outpt<-list(outres);names(outpt)<-"Results"
     return(outpt)
}
#dodo<-catchpertrip(intdir="C:/Temp",estdir="C:/Temp",species=8835250101,
#state=c(122),mode=c(3),wave=c(1:6),styr=1992,endyr=1992)





