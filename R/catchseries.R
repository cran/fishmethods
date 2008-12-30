catchseries<-function(estdir=NULL,species=NULL, state=NULL, byst=1, wave=NULL, bywave = 0 ,
mode=NULL, bymode=0, area=NULL, byarea=0, styr=NULL, endyr=NULL){
    if(is.null(estdir)) stop("Need main directory location of catch and effort files.")
    if(is.null(species)) stop("Need NODC code for species.")
    if(is.null(state)) {warning("No state code was specified. Data from all states will be used."); state<-c(1:56,72,78)}
    if(is.null(mode)) {warning("No mode code was specified. Data for all modes will be used.");mode<-c(1,2,3,4,5,6,7,8,9)}
    if(is.null(wave)) {warning("No wave code was specified. Data for all wave will be used.");wave<-c(1,2,3,4,5,6)}
    if(is.null(area)) {warning("No area code was specified. Data for all areas will be used.");area<-c(1:6)} 
   if(is.null(styr)) stop("Starting year is missing.")
    if(is.null(endyr)) stop("Ending year is missing.")

   
     if(length(grep("/",estdir))==1){
        dest<-ifelse(substr(estdir,nchar(estdir),nchar(estdir)) %in% c("/"),
          c(paste(estdir,"est",sep="")),c(paste(estdir,"/est",sep="")))
     }
    if(length(grep("\\\\",estdir))==1){
        dest<-ifelse(substr(estdir,nchar(estdir),nchar(estdir)) %in% c("\\"),
          c(paste(estdir,"est",sep="")),c(paste(estdir,"\\est",sep="")))
     }
if(styr>=1982 & endyr<=2004){
        if(!all(mode %in% c(3,6,7))) stop("Mode not valid for years selected.")
        }
if(styr>=2005){
        if(!all(mode %in% c(3,4,5,7))) stop("Mode not valid for years selected.")
        }
if(styr>=1982 & endyr>2004){
        if(!all(mode %in% c(3,4,5,6,7))) stop("Mode not valid for years selected.")
        }
if(any(state==121) & any(state==122)) stop("Use state 12 to combine Florida coasts.")
if(any(state==12) & any(state %in% c(121,122))) stop("You are mixing codes for Florida.")
stcnt<-0;spcnt<-0;modecnt<-0;areacnt<-0
outresults<-NULL
  species<-as.character(species)
	flag<-0
	for(yr in styr:endyr){
 	   for (j in 1:as.numeric(length(wave))){ 
      	  wv<-wave[j]     
         	ests<-read.csv(paste(dest,yr,"/","AG_",yr,wv,".csv",sep="")) 
          if(length(ests$SUB_REG)>0){
            flag<-flag+1
         	ests$SPCODE<-as.character(ests$SP_CODE)
            if(any(state %in% c(121,122))){
                      d3<-ests[ests$ST!=12,]
                      if(any(state==121)) d4<-ests[ests$ST==12 & ests$SUB_REG!=7,]
                      if(any(state==122)) d4<-ests[ests$ST==12 & ests$SUB_REG!=6,]
                      ests<-as.data.frame(rbind(d3,d4))
                      state1<-state[state!=121|state!=122]
                      state1<-c(state1,12)
                  }
            if(!any(state %in% c(121,122))) state1<-state
            if(any(ests$SP_CODE==species)) spcnt<-spcnt+1
            if(any(ests$ST %in% c(state1))) stcnt<-stcnt+1
            if(any(ests$MODE_FX %in% c(mode))) modecnt<-modecnt+1
            if(any(ests$AREA_X %in% c(area))) areacnt<-areacnt+1
                
	   	ests<-ests[ests$SP_CODE==species & ests$ST %in% c(state1) & ests$MODE_FX %in% c(mode) &
                      ests$AREA_X %in% c(area),]
            ests$TYPEA<-ests$ESTCLAIM
            ests$TYPEA.VAR<-ests$ESTCLVAR
            ests$TYPEB1<-ests$ESTHARV
            ests$TYPEB1.VAR<-ests$ESTHVAR
            ests$HARVEST<-ests$TYPEA+ests$TYPEB1
            ests$HARV.VAR<-ests$ESTHVAR+ests$ESTCLVAR
            ests$HARWGT<-ests$WGT_AB1
            ests$HARWGT.VAR<-ests$VAR_WAB1
            ests$RELEASED<-ests$ESTREL
            ests$REL.VAR<-ests$ESTRLVAR
            ests$TOTAL<-ests$HARVEST+ests$RELEASED
            ests$TOT.VAR<-ests$ESTCLVAR+ests$ESTHVAR+ests$ESTRLVAR

            ests<-subset(ests,select=c(ST,YEAR,WAVE,MODE_FX,AREA_X,SUB_REG,
                     TYPEA,TYPEA.VAR,TYPEB1,TYPEB1.VAR,HARVEST,HARV.VAR,RELEASED,
                     REL.VAR,HARWGT,HARWGT.VAR,TOTAL,TOT.VAR))
            if(flag==1) outresults<-ests
       	if(flag>1) outresults<-rbind(outresults,ests)
        }
     }
   }  
if(spcnt==0) stop("Species not found")       
if(stcnt==0) stop("State not found.")
if(areacnt==0) stop("Area not found.")
if(modecnt==0) stop("Mode not found.")
 
 if(flag>0){
   outresults<-outresults[!is.na(outresults$ST),]
   grps<-NULL;labels<-NULL
   if(byst==0){ grps[[1]]<-outresults[,2]
                names(grps)[1]<-"YEAR"
              }
   if(byst==1){ grps<-list(outresults[,1],outresults[,2])
                names(grps)[1:2]<-c("STATE","YEAR")
              }
   if(bywave==1){ grps[[length(grps)+1]]<-outresults[,3]
                  names(grps)[length(grps)]<-"WAVE"
                }
   if(bymode==1){ grps[[length(grps)+1]]<-outresults[,4]
 			 names(grps)[length(grps)]<-"MODE"
                }
   if(byarea==1){ grps[[length(grps)+1]]<-outresults[,5]
                  names(grps)[length(grps)]<-"AREA"
                }
   h1<-aggregate(outresults[,7:18],by=grps,sum,na.rm=T)
   lcnt<-length(grps)+1
   ccol<-ncol(h1)
   names(h1)[lcnt:ccol]<-c("TYPEA","TYPEA.VAR","TYPEB1","TYPEB1.VAR",
                 "HARVEST","HARV.VAR","RELEASED","REL.VAR",
                  "HARWGT_KG","HARWGT.VAR","TOTAL","TOT.VAR")
             h1$TYPEA.SE<-sqrt(h1$TYPEA.VAR)
             h1$TYPEB1.SE<-sqrt(h1$TYPEB1.VAR)
             h1$HARV.SE<-sqrt(h1$HARV.VAR)
             h1$REL.SE<-sqrt(h1$REL.VAR)
             h1$HARWGT.SE<-sqrt(h1$HARWGT.VAR)
             h1$TOT.SE<-sqrt(h1$TOT.VAR)   
             h1$TYPEA.PSE<-round(h1$TYPEA.SE/h1$TYPEA*100,1)
             h1$TYPEB1.PSE<-round(h1$TYPEB1.SE/h1$TYPEB1*100,1)
             h1$HARV.PSE<-round(h1$HARV.SE/h1$HARVEST*100,1)
             h1$REL.PSE<-round(h1$REL.SE/h1$RELEASED*100,1)
             h1$HARVWGT.PSE<-round(h1$HARWGT.SE/h1$HARWGT_KG*100,1)
             h1$TOT.PSE<-round(h1$TOT.SE/h1$TOTAL*100,1)
             h1<-h1[,!names(h1) %in% c("TYPEA.VAR","TYPEB1.VAR","HARV.VAR","REL.VAR","HARWGT.VAR","TOT.VAR")]
             details<-NULL
             details[1]<-c(paste("Species: ",species,sep=""))
             details[2]<-c(paste("State: ",paste(state1,collapse=","),sep=""))
             details[3]<-c(paste("Start and End Years: ",styr,"-",endyr,sep=""))
             details[4]<-c(paste("Wave: ",paste(wave,collapse=","),sep=""))
             details[5]<-c(paste("Mode: ",paste(mode,collapse=","),sep=""))
             details[6]<-c(paste("Area: ",paste(area,collapse=","),sep=""))
             outpt<-list(details,h1);names(outpt)<-c("Details","Results")
            return(outpt)
           }
          if(flag==0) return("No Data")
}
#dodo<-catchseries(estdir="C:/Temp",species=8835250101,
#state=c(12,12),byst=0,mode=c(3,4,5,7),bymode=0,wave=c(1:6),bywave=0,area=c(1:7),byarea=0, styr=2007,endyr=2007)




