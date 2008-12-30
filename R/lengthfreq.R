lengthfreq<-function(intdir=NULL, estdir=NULL,species=NULL, state=NULL, wave=NULL, mode=NULL, 
            area=NULL, styr=NULL, endyr=NULL, conveq=FALSE, parms=c(0,1)){
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
 if(!all(area %in% c(1,2,3,4,5,6,7))) stop ("Area not valid.")
 
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
if(!all(state %in% c(1,2,4,5,6,8,9,10,11,12,13,15:42,44:51,53:56))) stop ("Invalid state code.")

outresults<-NULL
  species<-as.character(species)
	flag<-0
	for(yr in styr:endyr){
 	   for (j in 1:as.numeric(length(wave))){ 
      	  wv<-wave[j] 
              t3<-read.csv(paste(din,yr,"/","I3_",yr,wv,".csv",sep=""))
              t3$ID_CODE<-as.character(t3$ID_CODE) 
              t3$SP_CODE<-as.character(t3$SP_CODE) 
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
              t3$MODE_FX %in% c(mode) & t3$AREA_X %in% c(area) & !is.na(t3$LNGTH),]
          
              t3<-subset(t3,select=c(ID_CODE,LNGTH))	

       	  if(length(t3$ID_CODE)==0){
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

             if(length(t1)==0) warning("State, mode or area not found." )
   
       	  m1<-merge(t1,t3,by.x="ID_CODE",by.y="ID_CODE",sort=T,all.y=T)  
  
       if(!is.na(m1$ID_CODE[1])|length(m1$ID_CODE)>1){
             flag<-flag+1
             if(conveq==TRUE) m1$LNGTH<-parms[1]+parms[2]*m1$LNGTH
             m1$INCH<-round(m1$LNGTH/25.4,3)
             m1$LCAT<-trunc(m1$INCH)
             lf<-as.data.frame(table(m1$ST,m1$MODE_FX,m1$AREA_X,m1$SUB_REG,m1$LCAT))
             names(lf)<-c("ST","MODE_FX","AREA_X","SUB_REG","LNGTH","NUMLEN")
             tr<-as.data.frame(table(m1$ST,m1$MODE_FX,m1$AREA_X,m1$SUB_REG))
             names(tr)<-c("ST","MODE_FX","AREA_X","SUB_REG","TOTAL")
             lens<-merge(lf,tr,by.x=c("ST","MODE_FX","AREA_X","SUB_REG"),
			by.y=c("ST","MODE_FX","AREA_X","SUB_REG"),sort=T)
             lens$RELFREQ<-ifelse(is.nan(lens$NUMLEN/lens$TOTAL),0,lens$NUMLEN/lens$TOTAL)
        
         	ests<-read.csv(paste(dest,yr,"/","AG_",yr,wv,".csv",sep="")) 
         	ests$SPCODE<-as.character(ests$SP_CODE)
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

	   	ests<-ests[ests$SP_CODE==species & ests$ST %in% c(state1) & ests$MODE_FX %in% c(mode) &
                      ests$AREA_X %in% c(area),]
            ests$ESTLAND<-ests$ESTCLAIM+ests$ESTHARV
            ests<-subset(ests,select=c(ST,WAVE,MODE_FX,AREA_X,SUB_REG,ESTLAND))
            totest<-sum(ests$ESTLAND)
            ffd<-merge(ests,lens,by.x=c("ST","MODE_FX","AREA_X","SUB_REG"),by.y=c("ST","MODE_FX","AREA_X","SUB_REG"),all.x=T,all.y=T)
            ffd<-ffd[!is.na(ffd$LNGTH),]
            ffd$NUMBER<-ifelse(is.nan(ffd$RELFREQ*ffd$ESTLAND)|is.na(ffd$RELFREQ*ffd$ESTLAND),0,ffd$RELFREQ*ffd$ESTLAND)
            if(flag==1) outresults<-ffd
       	if(flag>1) outresults<-rbind(outresults,ffd)
        }
   }
}         
  if(is.null(outresults)) stop ("No data found.")
            if(flag>0){
		f1<-aggregate(outresults$NUMBER,list(outresults$LNGTH),sum)
            f1$x<-(f1$x/sum(f1$x))*totest
            f1$PERCENT<-round(f1$x/sum(f1$x)*100,2)
		names(f1)<-c("INCH.GROUP","NUMBER","PERCENT") 
            d1<-as.data.frame(as.factor(seq(min(as.numeric(as.character(f1$INCH.GROUP))),
                  max(as.numeric(as.character(f1$INCH.GROUP))),1)))
            names(d1)<-"INCH.GROUP"
            f1<-merge(f1,d1,by.x="INCH.GROUP",by.y="INCH.GROUP",all.x=T,all.y=T)
            f1$NUMBER<-ifelse(is.na(f1$NUMBER),0,f1$NUMBER)
            f1$PERCENT<-ifelse(is.na(f1$PERCENT),0,f1$PERCENT)
            f1<-f1[order(as.numeric(as.character(f1$INCH.GROUP))),]
            details<-NULL
            details[1]<-c(paste("Species: ",species,sep=""))
            details[2]<-c(paste("State: ",paste(state,collapse=","),sep=""))
            details[3]<-c(paste("Start and End Years: ",styr,"-",endyr,sep=""))
            details[4]<-c(paste("Wave: ",paste(wave,collapse=","),sep=""))
            details[5]<-c(paste("Mode: ",paste(mode,collapse=","),sep=""))
            details[6]<-c(paste("Area: ",paste(area,collapse=","),sep=""))
            outpt<-list(details,f1);names(outpt)<-c("Details","Results")
            return(outpt)
           }
          if(flag==0) return("No Data")
}

#dodo<-lengthfreq(intdir="C:/Temp",estdir="C:/Temp",species=8835430101,
#state=c(25),mode=c(5),wave=c(4),area=c(1,2,3,4,5,7), styr=2007,
#endyr=2007,conveq=FALSE,parms=c(0,2))



