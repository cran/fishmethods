growth_sel<-function(age=NULL,size=NULL,weights=NULL,minlimit=NULL,maxlimit=NULL,minmax=NULL,
    switch_varpar=1,
    Linf=list(init=1000,lb=100,ub=2000,prior.mean=1000,prior.var=-0.5,prior.pdf=1),
    K=list(init=0.3,lb=0.1,ub=0.9,prior.mean=0.3,prior.var=-0.05,prior.pdf=1),
    t0=list(init=-0.5,lb=-2,ub=-0.0001,prior.mean=-0.5,prior.var=-0.5,prior.pdf=1),
    varpar=list(init=50.0,lb=10.0,ub=100.0,prior.mean=5,prior.var=-1.0,prior.pdf=1),
    tmb.control=list(maxit=5000,trace=F),nlminb.control=list(eval.max=100000,iter.max=1000),
    species_info=list(species=NULL,size_units=NULL)
    ){
  #Ensure all vectors are numeric
  types<-c(typeof(age),typeof(size),typeof(weights),typeof(minlimit),typeof(maxlimit),
      typeof(minmax))
  index<-which(types %in% c("character","logical"))
  if(length(index)>0) stop(paste("Non-numeric vectors: ",paste(c("age","size","weights","minlimit","maxlimit","minmax")[index],sep="",collapse=","),sep=""))
  
  #Check if minlimit and maxlimit are vectors or single values;if single, expand to age vector length
  if(length(minlimit)==1) minlimit<-rep(minlimit,length(age)) 
  if(length(maxlimit)==1) maxlimit<-rep(maxlimit,length(age)) 
  
  #Check if age, size, weights, minlimit, maxlimit, minmax vectors are the same length
   n_rows<-lengths(list(age,size,weights,minlimit,maxlimit,minmax))
   if(length(unique(n_rows))!=1) stop(paste("Vectors lengths not equal - age:",n_rows[1],
   ", size:",n_rows[2],", weights:",n_rows[3],", minlimit:",n_rows[4],", maxlimit:", n_rows[5],
   ", minmax:",n_rows[6],sep=""))
  
  #check switch_varpar only ranges 1-3
  if(!switch_varpar %in% c(1:3)) stop("switch_varpar must be set to 1,2 or 3")
  #check init is within bounds
  if(Linf$init<Linf$lb||Linf$init>Linf$ub) stop("Linf initial guess not inside boundary values (lb,ub)")
  if(K$init<K$lb||K$init>K$ub) stop("K initial guess not inside boundary values (lb,ub)")
  if(t0$init<t0$lb||t0$init>t0$ub) stop("t0 initial guess not inside boundary values (lb,ub)")
  if(varpar$init<varpar$lb||varpar$init>varpar$ub) stop("varpar initial guess not inside boundary values (lb,ub)")
  
  #check if parameter prior pdfs code in 1:4
  if(!Linf$prior.pdf %in% c(1:4)) stop("Linf prior.pdf value not valid")
  if(!K$prior.pdf %in% c(1:4)) stop("K prior.pdf value not valid")
  if(!t0$prior.pdf %in% c(1:4)) stop("t0 prior.pdf value not valid")
  if(!varpar$prior.pdf %in% c(1:4)) stop("varpar prior.pdf value not valid")
  
  #check pdf and priors are valid
  #lognormal
  if(Linf$prior.pdf==2 & Linf$prior.mean<=0.0) stop("Linf: Don't use a lognormal distribution for a negative prior")
  if(K$prior.pdf==2 & K$prior.mean<=0.0) stop("K: Don't use a lognormal distribution for a negative prior")
  if(t0$prior.pdf==2 & t0$prior.mean<=0.0) stop("t0: Don't use a lognormal distribution for a negative prior")
  if(varpar$prior.pdf==2 & varpar$prior.mean<=0.0) stop("varpar: Don't use a lognormal distribution for a negative prior")
  
   #Beta
  if(Linf$prior.pdf==4 & (Linf$prior.mean<=0.0||Linf$prior.mean>=1.0)) stop("Linf: Don't use a beta distribution for a prior outside (0,1)")
  if(K$prior.pdf==4 & (K$prior.mean<=0.0||K$prior.mean>=1.0)) stop("K: Don't use a beta distribution for a prior outside (0,1)")
  if(t0$prior.pdf==4 & (t0$prior.mean<=0.0||t0$prior.mean>=1.0)) stop("t0: Don't use a beta distribution for a prior outside (0,1)")
  if(varpar$prior.pdf==4 & (varpar$prior.mean<=0.0||varpar$prior.mean>=1.0)) stop("varpar: Don't use a beta distribution for a prior outside (0,1)")
  
  data.gs = list(
    age=age,    
    size = size,     
    weights=weights,  
    minlimit = minlimit,
    maxlimit = maxlimit,
    minmax = minmax,
    switch_varpar=switch_varpar,
    linf_prior=c(Linf$prior.mean,Linf$prior.var),
    linf_pdf=as.integer(Linf$prior.pdf),
    k_prior=c(K$prior.mean,K$prior.var),
    k_pdf=as.integer(K$prior.pdf),
    t0_prior=c(t0$prior.mean,t0$prior.var),
    t0_pdf=as.integer(t0$prior.pdf),
    varpar_prior=c(varpar$prior.mean,varpar$prior.var),
    varpar_pdf=as.integer(varpar$prior.pdf)
  )
   
  parameters.gs <- list( 
    linf = Linf$init,
    k = K$init,
    t0 = t0$init,
    varpar=varpar$init
  )     
  
  obj <- MakeADFun(data=c(model="grow_sel",data.gs),parameters.gs,DLL="fishmethods_TMBExports",inner.control=tmb.control,silent=TRUE) 
  opt<-nlminb(obj$par, obj$fn, obj$gr,lower=c(Linf$lb,K$lb,t0$lb,varpar$lb),
                 upper=c(Linf$ub,K$ub,t0$ub,varpar$ub),
                 control=nlminb.control)
  
  report_values<-obj$report()
  info_outpt<-data.frame(date=Sys.time(),title="VBGF Truncated Fit",
                         species=ifelse(is.null(species_info$species),"Unspecified",species_info$species),
                         model="VBGF",
                         size.units=ifelse(is.null(species_info$size_units),"Unspecified",species_info$size_units))
  
  likelihood_info<-data.frame(lk.fit=report_values$f_fit,lk.prior=report_values$f_prior,
                              lik.total=report_values$fval)
  convergence_info<-data.frame(convergence=opt$convergence,iterations=opt$iterations,
                               function_evaluations=as.numeric(opt$evaluations[1]),
                               gradient_evaluations=as.numeric(opt$evaluations[2]),
                               message=opt$message)
  
  parm_estimates<-as.data.frame(summary(sdreport(obj)))
  row.names(parm_estimates)[1:2]<-c("Linf","K")
  parm_estimates$Lower_Bound<-c(Linf$lb,K$lb,t0$lb,varpar$lb)
  parm_estimates$Upper_Bound<-c(Linf$ub,K$ub,t0$ub,varpar$ub)
  
  if(switch_varpar==1) row.names(parm_estimates)[4]<-"Sigma"
  if(switch_varpar==2) row.names(parm_estimates)[4]<-"CV"
  if(switch_varpar==3) row.names(parm_estimates)[4]<-"Variance to mean ratio"
  
  if(report_values$switch_filter!=1) message<-""
  if(report_values$switch_filter==1) message<-"One or more data entries were filtered because observed size < limit."
  data_outpt<-data.frame(age=age,obs_size=size,predicted_size=report_values$size_pred)
  outpt<-list(run_info=info_outpt,message=message,convergence_info=convergence_info,
                estimates=parm_estimates,likelihood=likelihood_info,predicted=data_outpt)
  return(outpt)
}
