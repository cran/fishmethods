pstrat<-function(intdir=NULL, estdir=NULL,pstdir=NULL, state=NULL,year=NULL, stwave=NULL, 
    endwave=NULL,psfactor=NULL){
    if(is.null(intdir)) stop("Need main directory location of intercept files.")
    if(is.null(estdir)) stop("Need main directory location of catch and effort files.")
    if(is.null(pstdir)) stop("Need main directory location to store poststratification files.")
    if(is.null(stwave)) warning("No starting wave code was specified.") 
    if(is.null(endwave)) warning("No ending wave code was specified.") 
    if(is.null(year))  stop("Year code is missing.")
   
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

      if(length(grep("/",pstdir))==1){
        pest<-ifelse(substr(pstdir,nchar(pstdir),nchar(pstdir)) %in% c("/"),
          c(paste(pstdir,"psest",sep="")),c(paste(estdir,"/psest",sep="")))
     }
    if(length(grep("\\\\",estdir))==1){
        pest<-ifelse(substr(pstdir,nchar(pstdir),nchar(pstdir)) %in% c("\\"),
          c(paste(pstdir,"psest",sep="")),c(paste(pstdir,"\\psest",sep="")))
     }
  dir.create(paste(pest,year,sep=""),showWarnings=F)
        sumrate<-NULL
        alltrips<-NULL
        numtype4<-NULL
        all2<-NULL
        all3<-NULL
        incadd<-NULL 
if(!all(state %in% c(1,2,4,5,6,8,9,10,11,12,13,15:42,44:51,53:56,90))) stop ("Invalid state code.")
#### Code equivalent to trpstrat.sas
       trp<-read.csv(paste(dest,year,"/","AG",year,".csv",sep="")) 
     trp<-trp[trp$POOL_FLG==1 & trp$OUTFLG==1 & trp$EX_FLG==1,] 
     trp$WAVE<-ifelse(trp$WAVE>=1 & trp$WAVE<2,1,
                  ifelse(trp$WAVE>=2 & trp$WAVE<3,2,
                  ifelse(trp$WAVE>=3 & trp$WAVE<4,3,
                  ifelse(trp$WAVE>=4 & trp$WAVE<5,4,
                  ifelse(trp$WAVE>=5 & trp$WAVE<6,5,6)))))
     trp<-trp[trp$WAVE %in% c(stwave:endwave) & trp$ST %in% c(state),]
     namep<-names(psfactor)[1:as.numeric(ncol(psfactor)-1)]
     
     
  for (j in stwave:endwave){      
# Code equivalent to trpstrat.sas, stats_ps.sas, and jonah.sas
       t1<-read.csv(paste(din,year,"/","I1_",year,j,".csv",sep=""))
       t1<-t1[t1$ST %in% c(state),]
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
       t1$ID_CODE<-as.character(t1$ID_CODE) 
       t1$MONTH<-as.numeric(substr(t1$ID_CODE,start=10,stop=11))
       t1$INT_NUM<-as.numeric(substr(t1$ID_CODE,start=14,stop=16))
       t1$YEAR<-as.numeric(substr(t1$ID_CODE,start=6,stop=9))
       if(nrow(t1)==0) next
       #REDUCE
        
        t1<-merge(t1,psfactor,by.x=namep,by.y=namep,all.x=T,all.y=T)
        t1$AREA.G<-as.character(t1$AREA.G)
        t1$AREA.G<-ifelse(is.na(t1$AREA.G),"XXX",t1$AREA.G)
        t1$AREA_S<-c(paste(t1$AREA.G,t1$AREA_X))
        t1<-subset(t1,select=c(YEAR, WAVE, SUB_REG, ST,MODE_FX,AREA_S,AREA_X,AREA,
              NUM_TYP4,HRSF,ADD_HRS,ID_CODE,CNTRBTRS,MONTH,INT_NUM,DIST))  

       if(year==1992 & j==1){
          t1$WAVE<-ifelse(t1$SUB_REG>=6 & t1$MONTH==1,1.1,
	          ifelse(t1$SUB_REG>=6 & t1$MONTH==2,1.2,t1$WAVE))
          t1$del<-ifelse(t1$SUB_REG>=6 & t1$MONTH>2,1,0)
          t1<-t1[t1$del==0,]
        }
       if(year==1988){
          if(j==4){
             t1$WAVE<-ifelse(t1$SUB_REG>=6 &t1$MONTH==7,4.1,
                    ifelse(t1$SUB_REG>=6 & t1$MONTH==8,4.2,t1$WAVE))
             t1$del<-ifelse(t1$SUB_REG>=6 & (t1$MONTH<7|t1$MONTH>8),1,0)
             t1<-t1[t1$del==0,]
           }
          if(j==5){
             t1$WAVE<-ifelse(t1$SUB_REG>=6 &t1$MONTH==9,5.1,
                  ifelse(t1$SUB_REG>=6 & t1$MONTH==10,5.2,t1$WAVE))
             t1$del<-ifelse(t1$SUB_REG>=6 & (t1$MONTH<9|t1$MONTH>10),1,0)
             t1<-t1[t1$del==0,]
           }
          if(j==6){
             t1$WAVE<-ifelse(t1$SUB_REG>=6 &t1$MONTH==11,6.1,
                ifelse(t1$SUB_REG>=6 & t1$MONTH==12,6.2,t1$WAVE))
             t1$del<-ifelse(t1$SUB_REG>=6 & (t1$MONTH<11|t1$MONTH>12),1,0)
             t1<-t1[t1$del==0,]
           }
       }
    t1<-t1[t1$INT_NUM<50,-c(14,15)]
    t1$CNTRBTRS<-ifelse(t1$CNTRBTRS %in% c(0,88,99,NA) & t1$NUM_TYP4==0,1,
                       t1$CNTRBTRS)
    t1$CNTRBTRS<-ifelse(t1$CNTRBTRS %in% c(0,88,99,NA) & t1$NUM_TYP4!=0,0,
                t1$CNTRBTRS) 
#for tripstats.sas  
     t1a<-t1
     t1a$AREA_X<-ifelse(t1a$AREA !=1,5,ifelse(t1a$AREA==1 & t1a$DIST==1,1,
              ifelse(t1a$AREA==1 & t1a$DIST==2,2,t1a$AREA_X)))
     rate<-as.data.frame(table(t1a$WAVE,t1a$SUB_REG,t1a$ST,t1a$MODE_FX,t1a$AREA_S))
         names(rate)<-c("WAVE","SUB_REG","ST","MODE_FX","AREA_S","PTTRIPS")
     rate<-rate[rate$PTTRIPS>0,] 
     sumrate<-rbind(sumrate,rate)
     t1<-t1[,-c(14)]

# For stats_ps.sas
       incomp<-t1
       incomp$ADD_HRS<-ifelse(is.na(incomp$ADD_HRS),0,incomp$ADD_HRS)
       incomp$HSUM<-ifelse(!is.na(incomp$HRSF)|!is.na(incomp$ADD_HRS),incomp$HRSF+incomp$ADD_HRS,0)
       incomp$N<-ifelse(!is.na(incomp$HRSF)|!is.na(incomp$ADD_HRS),1,0)
       incomp$CTRIP<-ifelse(incomp$ADD_HRS==0,1,0)
       incomp$CHRS<-ifelse(incomp$CTRIP==1,incomp$HRSF,0)
       incomp$INCTRIP<-ifelse(incomp$ADD_HRS>0,1,0)
       inc<-aggregate(incomp[,c(14:18)],list(incomp$YEAR,incomp$WAVE,incomp$SUB_REG,incomp$ST,
               incomp$MODE_FX, incomp$AREA_S),sum,na.rm=T)
         names(inc)[1:6]<-c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S") 
       inc$CELLTYPE<-ifelse(inc$INCTRIP==0,1,ifelse(inc$CTRIP==0,3,2))
       inc$FISH_ADJ<-ifelse(inc$INCTRIP==0,1,ifelse(inc$CTRIP==0,inc$HSUM/inc$N,
            inc$CHRS/inc$CTRIP))
      inc$CELLTYPE<-ifelse(is.na(inc$FISH_ADJ),4,inc$CELLTYPE)  
      inc$FISH_ADJ<-ifelse(is.na(inc$FISH_ADJ),1,inc$FISH_ADJ)
      inc<-inc[order(inc$YEAR,inc$WAVE,inc$SUB_REG,inc$ST,
                   inc$MODE_FX, inc$AREA_S),-c(7:11)]
      incadd<-rbind(incadd,inc)
      numtype<-t1[t1$CNTRBTRS==0 & t1$NUM_TYP4==0,c(12)]
      numtype<-as.data.frame(cbind(numtype,rep(1,length(numtype))))
      names(numtype)<-c("ID_CODE","TYPE4")
      numtype$ID_CODE<-as.character(numtype$ID_CODE)
      numtype4<-rbind(numtype,numtype)
      rm(numtype)

# for jonah_ps.sas
       tbl<-as.data.frame(table(t1$YEAR,t1$WAVE,t1$SUB_REG,t1$ST,t1$MODE_FX,t1$AREA_S))
       tbl<-tbl[tbl$Freq>0,]
       names(tbl)<-c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S","INT_TRIP")
       all_trip<-aggregate(t1$CNTRBTRS,list(t1$YEAR,t1$WAVE,t1$SUB_REG,t1$ST,t1$MODE_FX,t1$AREA_S),sum)
       names(all_trip)<-c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S","SMP_TRIP")
      all_trip<-merge(all_trip,tbl,by.x=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),
             by.y=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),all.x=T,all.y=T)
       alltrips<-rbind(alltrips,all_trip)
     rm(inc,t1,incomp,tbl,all_trip)
   
# for stats_ps.sas
       t2<-read.csv(paste(din,year,"/","I2_",year,j,".csv",sep=""))
       t2<-t2[t2$ST %in% c(state),]
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
       t2$ID_CODE<-as.character(t2$ID_CODE)
       t2$SP_CODE<-as.character(format(t2$SP_CODE,sci=F))
       t2$MONTH<-as.numeric(substr(t2$ID_CODE,start=10,stop=11))
       t2$INT_NUM<-as.numeric(substr(t2$ID_CODE,start=14,stop=16))
       t2$YEAR<-as.numeric(substr(t2$ID_CODE,start=6,stop=9))
     #REDUCE CODE
        t2<-merge(t2,psfactor,by.x=namep,by.y=namep,all.x=T,all.y=T)
        t2$AREA.G<-as.character(t2$AREA.G)
        t2$AREA.G<-ifelse(is.na(t2$AREA.G),"XXX",t2$AREA.G)
        t2$AREA_S<-c(paste(t2$AREA.G,t2$AREA_X))
        t2<-subset(t2,select=c(YEAR, WAVE, SUB_REG, ST,MODE_FX,AREA_S,AREA_X,
              ID_CODE,SP_CODE, NUM_FISH,DISPO, HRSF, ADD_HRS,CNTRBTRS,MONTH,INT_NUM))  
       if(year==1992 & j==1){
          t2$WAVE<-ifelse(t2$SUB_REG>=6 & t2$MONTH==1,1.1,
	          ifelse(t2$SUB_REG>=6 & t2$MONTH==2,1.2,t2$WAVE))
          t2$del<-ifelse(t2$SUB_REG>=6 & t2$MONTH>2,1,0)
          t2<-t2[t2$del==0,]
        }
       if(year==1988){
          if(j==4){
             t2$WAVE<-ifelse(t2$SUB_REG>=6 &t2$MONTH==7,4.1,
                    ifelse(t2$SUB_REG>=6 & t2$MONTH==8,4.2,t2$WAVE))
             t2$del<-ifelse(t2$SUB_REG>=6 & (t2$MONTH<7|t2$MONTH>8),1,0)
             t2<-t2[t2$del==0,]
           }
          if(j==5){
             t2$WAVE<-ifelse(t2$SUB_REG>=6 &t2$MONTH==9,5.1,
                  ifelse(t2$SUB_REG>=6 & t2$MONTH==10,5.2,t2$WAVE))
             t2$del<-ifelse(t2$SUB_REG>=6 & (t2$MONTH<9|t2$MONTH>10),1,0)
             t2<-t2[t2$del==0,]
           }
          if(j==6){
             t2$WAVE<-ifelse(t2$SUB_REG>=6 &t2$MONTH==11,6.1,
                ifelse(t2$SUB_REG>=6 & t2$MONTH==12,6.2,t2$WAVE))
             t2$del<-ifelse(t2$SUB_REG>=6 & (t2$MONTH<11|t2$MONTH>12),1,0)
             t2<-t2[t2$del==0,]
           }
       }
      t2<-t2[t2$INT_NUM<50,-c(15,16)]
      t2$DISPO<-ifelse(t2$DISPO %in% c(1,2),1,3)
      t2<-t2[order(t2$ID_CODE,t2$SP_CODE,t2$DISPO),]
      t2$RELEASE<-ifelse(t2$DISPO==1,t2$NUM_FISH,0)
      t2$HARVEST<-ifelse(t2$DISPO==3,t2$NUM_FISH,0)
       all2s<-as.data.frame(rowsum(t2[,15:16],c(paste(format(t2$ID_CODE,sci=F),format(t2$SP_CODE,sci=F)))))
     	 all2s$ID_CODE<-substr(as.character(rownames(all2s)),start=1,stop=16)
    	 all2s$SP_CODE<-substr(rownames(all2s),start=18,stop=27)
       all2s$HARVEST<-ifelse(is.na(all2s$HARVEST),0,all2s$HARVEST)
       all2s$RELEASE<-ifelse(is.na(all2s$RELEASE),0,all2s$RELEASE)
       all2s<-merge(t2[!duplicated(c(paste(t2$ID_CODE,t2$SP_CODE))),-c(10,11,15,16)],all2s,by.x=c("ID_CODE","SP_CODE"),
          by.y=c("ID_CODE","SP_CODE"),all.y=T)
       all2<-rbind(all2,all2s)
       all2<-all2[!is.na(all2$WAVE),]
       rm(all2s,t2)

  #TYPE 3
   t3<-read.csv(paste(din,year,"/","I3_",year,j,".csv",sep=""))
       t3<-t3[t3$ST %in% c(state),]
       t3$ID_CODE<-as.character(t3$ID_CODE)
       t3$SP_CODE<-as.character(format(t3$SP_CODE,sci=F))
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
       t3$MONTH<-as.numeric(substr(t3$ID_CODE,start=10,stop=11))
       t3$INT_NUM<-as.numeric(substr(t3$ID_CODE,start=14,stop=16))
       t3$YEAR<-as.numeric(substr(t3$ID_CODE,start=6,stop=9))
      #REDUCE CODE
        t3<-merge(t3,psfactor,by.x=namep,by.y=namep,all.x=T,all.y=T)
        t3$AREA.G<-as.character(t3$AREA.G)
        t3$AREA.G<-ifelse(is.na(t3$AREA.G),"XXX",t3$AREA.G)
        t3$AREA_S<-c(paste(t3$AREA.G,t3$AREA_X))
        t3<-subset(t3,select=c(YEAR, WAVE, SUB_REG, ST,MODE_FX,AREA_S,AREA_X,
              ID_CODE,SP_CODE, FSHINSP,WGT,LNGTH, HRSF, ADD_HRS,CNTRBTRS,MONTH,INT_NUM))  

       if(year==1992 & j==1){
          t3$WAVE<-ifelse(t3$SUB_REG>=6 & t3$MONTH==1,1.1,
	          ifelse(t3$SUB_REG>=6 & t3$MONTH==2,1.2,t3$WAVE))
          t2$del<-ifelse(t3$SUB_REG>=6 & t3$MONTH>2,1,0)
          t3<-t3[t3$del==0,]
        }
       if(year==1988){
          if(j==4){
             t3$WAVE<-ifelse(t3$SUB_REG>=6 &t3$MONTH==7,4.1,
                    ifelse(t3$SUB_REG>=6 & t3$MONTH==8,4.2,t3$WAVE))
             t3$del<-ifelse(t3$SUB_REG>=6 & (t3$MONTH<7|t3$MONTH>8),1,0)
             t3<-t3[t3$del==0,]
           }
          if(j==5){
             t3$WAVE<-ifelse(t3$SUB_REG>=6 &t3$MONTH==9,5.1,
                  ifelse(t3$SUB_REG>=6 & t3$MONTH==10,5.2,t3$WAVE))
             t3$del<-ifelse(t3$SUB_REG>=6 & (t3$MONTH<9|t3$MONTH>10),1,0)
             t3<-t3[t3$del==0,]
           }
          if(j==6){
             t3$WAVE<-ifelse(t3$SUB_REG>=6 & t3$MONTH==11,6.1,
                ifelse(t3$SUB_REG>=6 & t3$MONTH==12,6.2,t3$WAVE))
             t3$del<-ifelse(t3$SUB_REG>=6 & (t3$MONTH<11|t3$MONTH>12),1,0)
             t3<-t3[t3$del==0,]
           }
       }
       t3<-t3[t3$INT_NUM<50,-c(16,17)]
       t3a<-t3
       t3a$LCNT<-ifelse(is.na(t3a$LNGTH),0,1)
       t3a$WCNT<-ifelse(is.na(t3a$WGT),0,1)
       t3a$LNGTH<-ifelse(is.na(t3a$LNGTH),0,t3a$LNGTH)
       t3a$WGT<-ifelse(is.na(t3a$WGT),0,t3a$WGT)
       t3a$WGT2<-t3a$WGT^2
       t3a$LNGTH2<-t3a$LNGTH^2
       t3a$label<-c(paste(t3a$YEAR,t3a$WAVE,t3a$SUB_REG,t3a$ST,t3a$MODE_FX,
                 t3a$AREA_S, format(t3a$ID_CODE,sci=F),format(t3a$SP_CODE,sci=F))) 
       all3s<-as.data.frame(rowsum(t3a[,c(11,12,16:19)],t3a$label))
       all3s$label<-as.character(rownames(all3s))
       all3s$AVEWGT<-ifelse(all3s$WCNT==0,0,all3s$WGT/all3s$WCNT)
       all3s$AVELNGTH<-ifelse(all3s$LCNT==0,0,all3s$LNGTH/all3s$LCNT)    
       t3$label<-c(paste(t3a$YEAR,t3a$WAVE,t3a$SUB_REG,t3a$ST,t3a$MODE_FX,
                 t3a$AREA_S, format(t3a$ID_CODE,sci=F),format(t3a$SP_CODE,sci=F))) 
       all3s<-merge(t3[!duplicated(t3$label),-c(11,12)],all3s,
                by.x=c("label"),by.y=c("label"),all.x=T,all.y=T)

       all3<-rbind(all3,all3s)
       all3<-all3[!is.na(all3$WAVE),]
       rm(all3s,t3,t3a)

} #if j

### After loop through waves

#for trpstat.sas
    pstrips<-merge(trp,sumrate,by.x=c("WAVE","SUB_REG","ST","MODE_FX"),
             by.y=c("WAVE","SUB_REG","ST","MODE_FX"))
    sumtrips<-aggregate(pstrips$PTTRIPS,list(pstrips$WAVE,
             pstrips$SUB_REG,pstrips$ST,pstrips$MODE_FX),sum,na.rm=T)
           names(sumtrips)<-c("WAVE","SUB_REG","ST","MODE_FX","TOTAL")
    pstrips<-merge(pstrips,sumtrips,by.x=c("WAVE","SUB_REG","ST","MODE_FX"),
             by.y=c("WAVE","SUB_REG","ST","MODE_FX"))
    pstrips$PROP<-ifelse(is.nan(pstrips$PTTRIPS/pstrips$TOTAL),0,
                   pstrips$PTTRIPS/pstrips$TOTAL)
    pstrips$ESTRIPS<-pstrips$NUMRTRIP*pstrips$PROP
    pstrips$NUMVAR<-ifelse(pstrips$PROP==0,0,
                  pstrips$TOTALVAR*(pstrips$PROP^2)+
                   ((pstrips$PROP*(1-pstrips$PROP)*pstrips$NUMRTRIP^2)/pstrips$I_PLUS)-
                    (pstrips$PROP*(1-pstrips$PROP)*pstrips$TOTALVAR/pstrips$I_PLUS))

# for stat_ps.sas
 both23<-merge(all3,all2,by.x=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S","ID_CODE","SP_CODE"),
           by.y=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S","ID_CODE","SP_CODE"),all.x=T,all.y=T)
      both23$RELEASE<-ifelse(is.na(both23$RELEASE),0,both23$RELEASE)
      both23$HARVEST<-ifelse(is.na(both23$HARVEST),0,both23$HARVEST)
      both23$AREA_X<-ifelse(is.na(both23$AREA_X.x),both23$AREA_X.y,both23$AREA_X.x)
      both23$HRSF<-ifelse(is.na(both23$HRSF.x) & !is.na(both23$HRSF.y),both23$HRSF.y,
             ifelse(is.na(both23$HRSF.x) & is.na(both23$HRSF.y),0,both23$HRSF.x))
      both23$ADD_HRS<-ifelse(is.na(both23$ADD_HRS.x) & !is.na(both23$ADD_HRS.y),both23$ADD_HRS.y,
             ifelse(is.na(both23$ADD_HRS.x) & is.na(both23$ADD_HRS.y),0,both23$ADD_HRS.x))
      both23$CNTRBTRS<-ifelse(is.na(both23$CNTRBTRS.x) & !is.na(both23$CNTRBTRS.y),both23$CNTRBTRS.y,
             ifelse(is.na(both23$CNTRBTRS.x) & is.na(both23$CNTRBTRS.y),0,both23$CNTRBTRS.x))
      both23<-both23[,-c(9,10,12:14,23:26)]
      both23$FSHINSP<-ifelse(is.na(both23$FSHINSP),0,both23$FSHINSP)
       for(k in 1:ncol(both23)){
               if(is.numeric(both23[,k])) both23[is.na(both23[,k]),k] <- 0
              } 
      both23<-merge(both23,incadd,by.x=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),
         by.y=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),
         all.x=T)
      names(both23)[9]<-"CLAIM"

      both23$ADD_HRS<-ifelse(is.na(both23$ADD_HRS),0,both23$ADD_HRS)
      both23$CLAIM<-ifelse(both23$ADD_HRS>0,both23$FISH_ADJ*(both23$CLAIM/both23$HRSF),
                both23$CLAIM)
      both23$RELEASE<-ifelse(both23$ADD_HRS>0,both23$FISH_ADJ*(both23$RELEASE/both23$HRSF),
               both23$RELEASE)
      both23$HARVEST<-ifelse(both23$ADD_HRS>0,both23$FISH_ADJ*(both23$HARVEST/both23$HRSF),
               both23$HARVEST)
    
  #Add numtype for last section of stat_ps
       if(any(!is.na(match(both23$ID_CODE,numtype4$ID_CODE)))){
       both23<-merge(both23,numtype4,by.x=c("ID_CODE"),by.y=c("ID_CODE"),all.x=T)
       both23$CNTRBTRS<-ifelse(both23$TYPE4==1 & both23$CNTRBTRS==0,1,both23$CNTRBTRS)
       }
       both23$CLSSQR<-ifelse(both23$CNTRBTRS>1,both23$CNTRBTRS*((both23$CLAIM/both23$CNTRBTRS)^2),
                      both23$CLAIM^2)    
       both23$RLSSQR<-both23$RELEASE^2
       both23$HLSSQR<-both23$HARVEST^2
   
   #Output data file for stat_ps.sas 
     data<-as.data.frame(rowsum(both23[,c(9:19,26:28)], 
          c(paste(format(both23$YEAR,width=max(nchar(both23$YEAR))),
                  format(both23$WAVE,width=max(nchar(both23$WAVE))),
                  format(both23$SUB_REG,width=max(nchar(both23$SUB_REG))),
                  format(both23$ST,width=max(nchar(both23$ST))),
                  format(both23$MODE_FX,width=max(nchar(both23$MODE_FX))),
                  format(both23$AREA_S,width=max(nchar(both23$AREA_S))),
                  format(both23$SP_CODE,sci=F,width=max(nchar(both23$SP_CODE)))))),na.rm=T)

      names(data)<-c("TSPCLAIM","TSP_WGT","TSP_LEN","TSP_LEX","TSP_EXAM",
                "WGTSSQ","LNGSSQ","TSPAVEW","AVELNGTH","TSP_REL","TSP_HARV",
                 "CLSSQR","RLSSQR","HLSSQR")
          
      data$TSPAVEL<-ifelse(data$TSP_LEX>0,data$TSP_LEN/data$TSP_LEX,0)
      data$TSPAVEW<-ifelse(data$TSPCLAIM==0,0,data$TSP_WGT/data$TSP_EXAM)
      data<-data[,-c(9)]
      data$label<-rownames(data)
       
      both23$label<- c(paste(format(both23$YEAR,width=max(nchar(both23$YEAR))),
                  format(both23$WAVE,width=max(nchar(both23$WAVE))),
                  format(both23$SUB_REG,width=max(nchar(both23$SUB_REG))),
                  format(both23$ST,width=max(nchar(both23$ST))),
                  format(both23$MODE_FX,width=max(nchar(both23$MODE_FX))),
                  format(both23$AREA_S,width=max(nchar(both23$AREA_S))),
                  format(both23$SP_CODE,sci=F,width=max(nchar(both23$SP_CODE)))))
      b<-unique(both23[,c(1:6,8,29)])
      data<-merge(data,b,by.x="label",by.y="label",all.x=T)
      rm(b)
      data<-data[order(data$YEAR,data$WAVE,data$SUB_REG,data$ST,data$MODE_FX,
                data$AREA_S,data$SP_CODE),-c(1)]
       wvs<-as.numeric(c(levels(as.factor(data$WAVE))))
         for(o in 1:as.numeric(length(wvs))){
             ww<-wvs[o]
             write.table(data[data$WAVE==ww,],file=paste(pest,year,"/","STS",year,ww,".csv",sep=""),sep=",",row.names=F)
         }


# For jonah.sas
        allocate<-pstrips[pstrips$POOL_FLG==1 & pstrips$OUTFLG==1,]
        allocate$AREA_S<-as.character(allocate$AREA_S)

        pre<-merge(data,alltrips,by.x=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),
               by.y=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),all.y=T,all.x=T)

        pre<-merge(pre,allocate,by.x=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),
               by.y=c("YEAR","WAVE","SUB_REG","ST","MODE_FX","AREA_S"),all.y=T,all.x=T)
         
        pre$F_PER_T<-pre$TSPCLAIM/pre$SMP_TRIP
        pre$UNARPTRP<-pre$TSP_REL/pre$INT_TRIP
        pre$UNAHPTRP<-pre$TSP_HARV/pre$INT_TRIP
        pre$ESTCLAIM<-pre$F_PER_T*pre$ESTRIPS
        pre$ESTREL<-pre$UNARPTRP*pre$ESTRIPS
        pre$ESTHARV<-pre$UNAHPTRP*pre$ESTRIPS
        pre$TOT_CAT<-pre$ESTCLAIM+pre$ESTREL+pre$ESTHARV
        pre$LANDING<-pre$ESTCLAIM+pre$ESTHARV
        pre$VARLNGTH<-ifelse(pre$TSP_LEX>1,(pre$LNGSSQ-((pre$TSP_LEN^2)/pre$TSP_LEX))/(pre$TSP_LEX*(pre$TSP_LEX-1)),0)
        pre$VARWGT<-ifelse(pre$TSP_EXAM>1,(pre$WGTSSQ-((pre$TSP_WGT^2)/pre$TSP_EXAM))/(pre$TSP_EXAM*(pre$TSP_EXAM-1)),0)
        pre$VARCLAIM<-ifelse(pre$SMP_TRIP>1,(pre$CLSSQR-((pre$TSPCLAIM^2)/pre$SMP_TRIP))/(pre$SMP_TRIP*(pre$SMP_TRIP-1)),0)
        pre$VAREL<-ifelse(pre$SMP_TRIP>1,(pre$RLSSQR-((pre$TSP_REL^2)/pre$INT_TRIP))/(pre$INT_TRIP*(pre$INT_TRIP-1)),0)
        pre$VARHARV<-ifelse(pre$SMP_TRIP>1,(pre$HLSSQR-((pre$TSP_HARV^2)/pre$INT_TRIP))/(pre$INT_TRIP*(pre$INT_TRIP-1)),0)
        pre$ESTCLVAR<-ifelse(pre$VARCLAIM>0,(pre$VARCLAIM*pre$ESTRIPS^2)+(pre$NUMVAR*pre$F_PER_T^2)-(pre$VARCLAIM*pre$NUMVAR),0)
        pre$ESTRLVAR<-ifelse(pre$VAREL>0,(pre$VAREL*pre$ESTRIPS^2)+(pre$NUMVAR*pre$UNARPTRP^2)-(pre$VAREL*pre$NUMVAR),0)
        pre$ESTHVAR<-ifelse(pre$VARHARV>0,(pre$VARHARV*pre$ESTRIPS^2)+(pre$NUMVAR*pre$UNAHPTRP^2)-(pre$VARHARV*pre$NUMVAR),0)
        pre$ESTWGT<-ifelse(pre$TSP_EXAM>0,(pre$TSP_WGT/pre$TSP_EXAM)*pre$ESTCLAIM,0)
        pre$ESTWTVAR<-ifelse(pre$TSP_EXAM>0,(pre$VARWGT*pre$ESTCLAIM^2)+(pre$VARCLAIM*(pre$TSPAVEW^2))-
                       (pre$VARCLAIM*pre$VARWGT),0)
        pre$WGT_AB1<-ifelse(pre$ESTCLAIM>0,pre$TSPAVEW*pre$LANDING,0)
        pre$LBS_AB1<-pre$WGT_AB1*2.2046226
        pre$VAR_WAB1<-ifelse(pre$ESTCLAIM>0,(pre$VARWGT*pre$LANDING^2)+(pre$TSPAVEW^2*(pre$ESTCLVAR+pre$ESTHVAR))-
                       (pre$VARWGT*(pre$ESTCLVAR+pre$ESTHVAR)),0)
        pre$VAR_LBS<-pre$VAR_WAB1*4.8603608
        pre$TOT_VAR<-pre$ESTCLVAR+pre$ESTHVAR+pre$ESTRLVAR
        pre$LAND_VAR<-pre$ESTCLVAR+pre$ESTHVAR
        pre<-pre[,-c(17:19,30:35,37:40,43:49)]
        x1<-aggregate(pre[,c(8,11,12)],pre[,c(1,2,3,4,18)],sum)
           names(x1)[6:8]<-c("ST_WGT","ST_EXAM","ST_SSQ")
           x1$ST_AVE<-ifelse(x1$ST_EXAM>0,x1$ST_WGT/x1$ST_EXAM,0)
           x1$STVAR<-ifelse(x1$ST_EXAM>1,((x1$ST_SSQ-((x1$ST_WGT^2)/x1$ST_EXAM))/(x1$ST_EXAM-1)),0)
           x1$FLAG1<-1
         x2<-aggregate(pre[,c(8,11,12)],pre[,c(1,2,3,18)],sum)
           names(x2)[5:7]<-c("SUB_WGT","SUB_EXAM","SUB_SSQ")
           x2$SUB_AVE<-ifelse(x2$SUB_EXAM>0,x2$SUB_WGT/x2$SUB_EXAM,0)
           x2$SUBVAR<-ifelse(x2$SUB_EXAM>1,((x2$SUB_SSQ-((x2$SUB_WGT^2)/x2$SUB_EXAM))/(x2$SUB_EXAM-1)),0)
           x2$FLAG2<-2
         pre<-merge(pre,x2,by.x=c("SP_CODE","YEAR","WAVE","SUB_REG"),
                        by.y=c("SP_CODE","YEAR","WAVE","SUB_REG"),all.x=T,all.y=T)
         pre<-merge(pre,x1,by.x=c("SP_CODE","YEAR","WAVE","SUB_REG","ST"),
                        by.y=c("SP_CODE","YEAR","WAVE","SUB_REG","ST"),all.x=T,all.y=T)
         pre$FLAG_WGT<-0
         pre$ESTCLAIM<-ifelse(is.na(pre$ESTCLAIM),0,pre$ESTCLAIM)
         pre$SUB_EXAM<-ifelse(is.na(pre$SUB_EXAM),0,pre$SUB_EXAM)
         pre$SUBVAR<-ifelse(is.na(pre$SUBVAR),0,pre$SUBVAR)
         pre$ST_EXAM<-ifelse(is.na(pre$ST_EXAM),0,pre$ST_EXAM)
         pre$STVAR<-ifelse(is.na(pre$STVAR),0,pre$STVAR)
         pre$VARWGT<-ifelse(is.na(pre$VARWGT),0,pre$VARWGT)
         pre$ESTHARV<-ifelse(is.na(pre$ESTHARV),0,pre$ESTHARV)
         pre$TSP_EXAM<-ifelse(is.na(pre$TSP_EXAM),0,pre$TSP_EXAM)
         pre$BOOL1<-ifelse(pre$ESTCLAIM==0 & pre$ESTHARV>0,1,0)
         pre$BOOL2<-ifelse(pre$ESTCLAIM>0 & pre$TSP_EXAM<=1,1,0)
         pre$BOOL2<-ifelse(pre$ESTCLAIM>0 & pre$VARWGT==0,1,pre$BOOL2)

   for(k in 1:as.numeric(length(pre$YEAR))){
             if(pre$BOOL1[k]==1 | pre$BOOL2[k]==1){
             if(pre$SUB_REG[k]<3){
                if(pre$SUB_EXAM[k]>1 & pre$SUBVAR[k]>0){
                    pre$TSP_EXAM[k]<-pre$SUB_EXAM[k]
                    pre$TSP_WGT[k]<-pre$SUB_WGT[k]
                    pre$VARWGT[k]<-pre$SUBVAR[k]
                    pre$FLAG_WGT[k]<-pre$FLAG2[k]
                 }
                if(pre$SUB_EXAM[k]<=1 & pre$SUBVAR[k]==0){
                    pre$TSP_WGT[k]<-pre$ST_WGT[k] 
                    pre$TSP_EXAM[k]<-pre$ST_EXAM[k] 
                    pre$VARWGT[k]<-pre$STVAR[k] 
                    pre$FLAG_WGT[k]<-pre$FLAG1[k]
                 }
            }
            if(pre$SUB_REG[k]>=3){
               if(pre$ST_EXAM[k]>1 & pre$STVAR[k]>0){
                    pre$TSP_EXAM[k]<-pre$ST_EXAM[k]
                    pre$TSP_WGT[k]<-pre$ST_WGT[k]
                    pre$VARWGT[k]<-pre$STVAR[k]
                    pre$FLAG_WGT[k]<-pre$FLAG1[k]
                }
               if(pre$ST_EXAM[k]<=1 & pre$STVAR[k]==0){
 			  pre$TSP_EXAM[k]<-pre$SUB_EXAM[k]
                    pre$TSP_WGT[k]<-pre$SUB_WGT[k]
                    pre$VARWGT[k]<-pre$SUBVAR[k]
                    pre$FLAG_WGT[k]<-pre$FLAG2[k]
                }
            }
        }
    }

   pre$FLAG_WGT<-ifelse(pre$TSP_EXAM<=1,9,pre$FLAG_WGT)
   pre$TSPAVEW<-ifelse(pre$FLAG_WGT %in% c(1,2),pre$TSP_WGT/pre$TSP_EXAM,pre$TSPAVEW)
   pre$ESTWGT<-ifelse(pre$FLAG_WGT %in% c(1,2),pre$TSPAVEW*pre$ESTCLAIM,pre$ESTWGT)
   pre$WGT_AB1<-ifelse(pre$FLAG_WGT %in% c(1,2),pre$TSPAVEW*pre$LANDING,pre$WGT_AB1)
   pre$LBS_AB1<-pre$WGT_AB1*2.2046226
   pre$ESTWTVAR<-ifelse(pre$FLAG_WGT %in% c(1,2),(pre$VARWGT*pre$ESTCLAIM^2)+
                   (pre$VARCLAIM*(pre$TSPAVEW^2))-(pre$VARCLAIM*pre$VARWGT),
                   pre$ESTWTVAR)
   pre$VAR_WAB1<-ifelse(pre$FLAG_WGT %in% c(1,2),(pre$VARWGT*(pre$LANDING^2))+
                   ((pre$TSPAVEW^2)*(pre$ESTCLVAR+pre$ESTHVAR))-
                    (pre$VARWGT*(pre$ESTCLVAR+pre$ESTHVAR)),
                     pre$VAR_WAB1)
   pre$VAR_LBS<-pre$VAR_WAB1*4.8603608
   pre$TTRIP<-pre$PTTRIPS
   pre<-pre[,-c(27:34,37:46,79,85,87:88)] 
   pre<-pre[order(pre$SP_CODE,pre$SUB_REG,pre$YEAR,pre$WAVE,pre$ST,pre$MODE_FX,pre$AREA_S),]
   pre$DATE1<-date()
      for(k in 1:ncol(pre)){
               if(is.numeric(pre[,k])) pre[is.na(pre[,k]),k] <- 0
              } 
#output as in SAS program, as separate wave files 
         wvs<-as.numeric(c(levels(as.factor(pre$WAVE))))
         for(o in 1:as.numeric(length(wvs))){
             ww<-wvs[o]
             write.table(pre[pre$WAVE==ww,],file=paste(pest,year,"/","AG_",year,ww,".csv",sep=""),sep=",",row.names=F)
         }

}


