##############################################################################
#
#    Program: extractMRFSS.R
#    
#    CREATE DIRECTORIES FOR EACH YEAR OF DATA, EXTRACT MRFSS INTERCEPT OR
#      CATCH/EFFORT MASTER FILES FROM SAS TRANSPORT (.xpt) FILES AND SAVE 
#      AS COMMA_DELIMITED FILES (.csv) IN SPECIFIED DIRECTORY 
#
#    Gary Nelson
#    Massachusetts Division of Marine Fisheries
#    30 Emerson Avenue
#    Gloucester, MA 01930
#    gary.nelson@state.ma.us
#
#############################################################################

extractMRFSS<-function(indir=NULL,outdir=NULL,type=NULL,state=NULL,styr=NULL,endyr=NULL){
    require(foreign)
    if(is.null(indir)) stop("Need directory location for .xpt files.")
    if(is.null(styr)) stop("Starting year is missing.")
    if(is.null(endyr)) stop("Ending year is missing.")
    if(is.null(type)) stop("Data type is missing.")
    if(is.null(outdir)) outdir<-indir

 if(type==1){
      if(length(grep("/",indir))==1){
        din<-ifelse(substr(indir,nchar(indir),nchar(indir)) %in% c("/"),
          c(paste(indir,"int",sep="")),c(paste(indir,"/int",sep="")))
       }
      if(length(grep("\\\\",indir))==1){
        din<-ifelse(substr(indir,nchar(indir),nchar(indir)) %in% c("\\"),
          c(paste(indir,"int",sep="")),c(paste(indir,"\\int",sep="")))
      }
     if(length(grep("/",outdir))==1){
        dout<-ifelse(substr(outdir,nchar(outdir),nchar(outdir)) %in% c("/"),
          c(paste(outdir,"int",sep="")),c(paste(outdir,"/int",sep="")))
      }
     if(length(grep("\\\\",outdir))==1){
       dout<-ifelse(substr(outdir,nchar(outdir),nchar(outdir)) %in% c("\\"),
          c(paste(outdir,"int",sep="")),c(paste(outdir,"\\int",sep="")))
      }
  
    for(i in styr:endyr){
        dir.create(paste(dout,i,sep=""))
        data<-read.xport(paste(din,i,"ag.xpt",sep=""))
 	  for(j in 1:as.numeric(length(data))){
            outpt<-as.data.frame(data[j][[1]])
            if(!is.null(state)) outpt<-outpt[outpt$ST %in% c(state),]
            write.table(outpt,file=paste(dout,i,"/",names(data[j]),".csv",
              sep=""),sep=",",row.names=F)
    	    }
      }
  }

 if(type==2){
    if(length(grep("/",indir))==1){
       din<-ifelse(substr(indir,nchar(indir),nchar(indir)) %in% c("/"),
       c(paste(indir,"est",sep="")),c(paste(indir,"/est",sep="")))
     }
    if(length(grep("\\\\",indir))==1){
       din<-ifelse(substr(indir,nchar(indir),nchar(indir)) %in% c("\\"),
       c(paste(indir,"est",sep="")),c(paste(indir,"\\est",sep="")))
     }
    if(length(grep("/",outdir))==1){
       dout<-ifelse(substr(outdir,nchar(outdir),nchar(outdir)) %in% c("/"),
       c(paste(outdir,"est",sep="")),c(paste(outdir,"/est",sep="")))
     }
    if(length(grep("\\\\",outdir))==1){
       dout<-ifelse(substr(outdir,nchar(outdir),nchar(outdir)) %in% c("\\"),
       c(paste(outdir,"est",sep="")),c(paste(outdir,"\\est",sep="")))
     }
    for(i in styr:endyr){
        dir.create(paste(dout,i,sep=""),showWarnings=F)
        data<-read.xport(paste(din,i,"ag.xpt",sep=""))
 	  for(j in 1:as.numeric(length(data))){
            outpt<-as.data.frame(data[j][[1]])
            if(!is.null(state)) outpt<-outpt[outpt$ST %in% c(state),]
            write.table(outpt,file=paste(dout,i,"/",names(data[j]),".csv",
             sep=""),sep=",",row.names=F)
    	    }
     }
 }

}
#extractMRFSS(indir="C:\\temp\\",outdir="C:/temp",type=2,state=25,styr=1982,endyr=1982)
