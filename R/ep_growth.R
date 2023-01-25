ep_growth<-function(len=NULL,age=NULL,Nh=NULL,nh=NULL,starts=list(Linf=60,k=0.1,a0=-0.01,CV=0.5),bin_size=2,
                    nlminb.control=list(eval.max=5000,iter.max=5000,trace=10),
                    tmb.control=list(maxit=5000,trace=FALSE),plot=TRUE){
 #Check for missing data
  if(any(is.na(len))) stop("Missing value present in len")
  if(any(is.na(age))) stop("Missing value present in age")
  if(any(is.na(Nh))) stop("Missing value present in Nh")
  if(any(is.na(nh))) stop("Missing value present in Nh")
  
 #Check for different vector lengths
  vl<-c(length(len),length(age),length(Nh),length(nh))
  if(!all(vl))  stop("Unequal vector lengths in input data")
  data_EP<-data.frame(length=len,age=age,Nh=Nh,nh=nh)
  
  #to store results
  results<-matrix(NA, nrow = 4, ncol = 1)  
  colnames(results)<-c("EP")
  rownames(results)<-c("Linf","k","CV","a0")

 #to define variables
  data_EP$bin<-((data_EP$length)%/%bin_size)*bin_size
  data_EP$rlength<-data_EP$bin
  bin_max<-max(data_EP$bin)*2 
  len_max<-bin_max + bin_size
  len_min<-0

###reshape data for estimation
#extract unique bin sizes and calculate EP weight
  tdat<-data_EP[!duplicated(data_EP$bin),]
  tdat<-tdat[order(tdat$bin),]
  tdat$rlength<-tdat$bin
  tdat$ep<-tdat$nh/tdat$Nh

#to get table of Nh/nh
  tdat_all<- data.frame(rlength=seq(len_min,bin_max,bin_size))
  tdat_all$Nh<-0
  tdat_all$nh<-0
  tdat_all$ep<-1
  for (i_all in 1:length(tdat_all$rlength)){
  if(is.element(tdat_all$rlength[i_all],tdat$rlength)){
    temp_all<-tdat[tdat$rlength==tdat_all$rlength[i_all],]
    tdat_all$Nh[i_all]<-temp_all$Nh
    tdat_all$nh[i_all]<-temp_all$nh
    tdat_all$ep[i_all]<-temp_all$ep
   }
 }

#to get counts of lengths in each bin
 sel.tab<-as.data.frame(table(data_EP$length,data_EP$age))
 names(sel.tab)<-c("length","age","n")
 sel.tab<-sel.tab[sel.tab$n!=0,]
 sel.tab$age<-as.numeric(as.character(sel.tab$age))
 sel.tab$length<-as.numeric(as.character(sel.tab$length))
 sel.tab$rlength<-((sel.tab$length)%/%bin_size)*bin_size

###hooray! ready for estimation!
#EP estimation
tmb.data.EP<-list(
  age=sel.tab$age,    
  vlength = sel.tab$length,     
  vn = sel.tab$n,  
  rlength = tdat_all$rlength,
  vp = tdat_all$ep,
  age_ind = 1:(max(sel.tab$age)),
  len_ind = len_min:len_max,
  sample_size = length(data_EP$age),
  bin_size = bin_size
)

#starting values for parameter estimates
parameters.EP <- list( 
  log_Linf = log(starts$Linf),
  log_k = log(starts$k),
  log_CV = log(starts$CV),
  ao = starts$a0,
  lambda = rep(0,length(tmb.data.EP$age_ind)-1)
)     


obj.EP <- TMB::MakeADFun(data=c(model="EP_likelihood",tmb.data.EP),parameters.EP, 
                                   DLL="fishmethods_TMBExports", 
                                   inner.control=tmb.control)

#optimizer
opt.EP<-nlminb(obj.EP$par, obj.EP$fn, obj.EP$gr, 
               control = nlminb.control) 

ad_obj_EP<-sdreport(obj.EP, getReportCovariance = T)

outparms<-summary(ad_obj_EP)
index<-which(rownames(outparms)=="Linf")
output_parms<-outparms[c(index:nrow(outparms)),]
outs<-matrix("NA",nrow=6,ncol=1)
dimnames(outs)<-list(c("objective","convergence","iterations","function_eval","function_gradient","message"))
outs[1,1]<-opt.EP$objective;outs[2,1]<-opt.EP$convergence;outs[3,1]<-opt.EP$iterations;outs[4,1]<-as.numeric(opt.EP$evaluations[1]);
outs[5,1]<-as.numeric(opt.EP$evaluations[2]);outs[6,1]<-opt.EP$message
outs<-data.frame(outs)
names(outs)<-"value"


###plot results
#VonB function
VonB<-function(Linf,K,ao,age){
  Length = Linf*(1-exp(-K*(age-ao)))
  return(Length)
}
age.all<-1:max(data_EP$age)
plot.EP<-VonB(exp(opt.EP$par[1]),exp(opt.EP$par[2]),opt.EP$par[4],age.all)
results<-list(model=outs,parameter_estimates=output_parms,predicted=VonB(exp(opt.EP$par[1]),exp(opt.EP$par[2]),opt.EP$par[4],data_EP$age))

if(plot){
    plot(length~age,data=data_EP,pch=16,col="red",ylab="Length",xlab="Age",cex=1.2,cex.lab=1.2,cex.axis=1.2)
    lines(plot.EP~age.all,lwd=1.7)
    legend("topleft",bty="n",legend=c("Observed","Predicted"),lwd=1.7,cex=1.2,lty=c(NA,1),col=c("red","black"),pch=c(16,NA))
}  
return(results)
}#function end


