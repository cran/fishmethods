growth_LEP<-function(l1=NULL,l2=NULL,dt=NULL,measurer = NULL,
                      gmodel=1,use_parameter_boundaries=T,graphs=T,
                      K_start_bounds=list(K1=NULL,lower_K1=0,upper_K1=Inf,K2=NULL,lower_K2=0,upper_K2=Inf), 
                      mu_Linf_start_bounds=list(mu_Linf=NULL,lower_mu_Linf=0,upper_mu_Linf=Inf, sigma_mu_Linf=NULL, lower_sigma_mu_Linf=0, 
                                upper_sigma_mu_Linf=Inf),
                      A_start_bounds=list(mean_age=NULL,lower_mean_age=0,upper_mean_age=Inf,sigma_age=NULL,lower_sigma_age=0,
                               upper_sigma_age=Inf),
                      resid_error_start_bounds=list(sigma_resid=NULL,lower_sigma_resid=0,upper_sigma_resid=Inf),
                      measurer_error_start_bounds=list(use_measurer=F,sigma_measure=NULL,lower_sigma_measure=0,upper_sigma_measure=Inf),
                      vb_log_k_parms=list(alpha=NULL,lower_alpha=0,upper_alpha=Inf,fix_beta=T,beta=NULL,lower_beta=0,upper_beta=Inf),
                      nlminb.control=list(eval.max=10000,iter.max=10000,trace=10),
                      tmb.control=list(maxit=10000,trace=FALSE)){

#Error checks
if(is.null(l1)|is.null(l2)|is.null(dt)) stop("NULL for l1, l2 or dt")
if(length(l1)!=length(l2)) stop("l1 and l2 are different lengths")
if(!is.vector(dt)) stop("dt is not a vector")
  
if(length(dt)!=length(l1)) stop("dt and l1 are different lengths")
if(measurer_error_start_bounds$use_measurer==T){
    if(is.null(measurer)) stop("NULL for measurer ")
    if(!is.vector(measurer)) stop("measurer is not a vector")
    if(!length(measurer) %in% c(length(l1),length(l2),length(dt))) stop("measurer different length than l1, l2 or dt")
    if(is.null(measurer_error_start_bounds$sigma_measure)) stop("Missing start value for sigma_measure.")
  }  
  
if(!gmodel %in% c(1,2)) stop("model number invalid")
if(gmodel==1 & is.null(K_start_bounds$K1)) stop("Missing start value for K1.")
if(gmodel==2){
    if(is.null(K_start_bounds$K1)|is.null(K_start_bounds$K2)) stop("Missing start value for K1 or K2.")
    if(is.null(vb_log_k_parms$alpha)) stop("Missing start value for alpha.")
    if(is.null(vb_log_k_parms$beta)) stop("Missing start value for beta.")
}
if(is.null(mu_Linf_start_bounds$mu_Linf)) stop("Missing start value for mu.")
if(is.null(mu_Linf_start_bounds$sigma_mu_Linf)) stop("Missing start value for sigma_mu.")
if(is.null(A_start_bounds$mean_age)) stop("Missing start value for mean_age.")
if(is.null(A_start_bounds$sigma_age)) stop("Missing startvalue for sigma_age.")
if(is.null(resid_error_start_bounds$sigma_resid)) stop("Missing start value for sigma_resid.")
if(measurer_error_start_bounds$use_measurer==F) measurer<-rep(0,length(l1))

## TMB data and parameters
map<-list()
data <- list(l1=l1,
             l2=l2,
             dt=dt,
             measure=measurer,
             model1=gmodel,
             fix_beta=sum(vb_log_k_parms$fix_beta),
             model_measure=sum(measurer_error_start_bounds$use_measurer))


if(gmodel==1){
   map<-list(log_K2=factor(NA),log_alpha=factor(NA),log_beta=factor(NA))
   if(measurer_error_start_bounds$use_measurer==F){
     map<-append(map,list(log_measure=factor(NA)))
     mm<-1e-15
   }
   if(measurer_error_start_bounds$use_measurer==T){
     mm<-measurer_error_start_bounds$sigma_measure
   }
   
   parameters<-list(
    log_K1=log(K_start_bounds$K1),
    log_K2=log(1e-15),
    log_mu_Linf=log(mu_Linf_start_bounds$mu_Linf),
    log_mean_age=log(A_start_bounds$mean_age),
    log_sigma_mu_Linf=log(mu_Linf_start_bounds$sigma_mu_Linf),
    log_sigma_age=log(A_start_bounds$sigma_age),
    log_sigma_resid=log(resid_error_start_bounds$sigma_resid),
     A=rep(A_start_bounds$mean_age,length(l1)),#age random effects
    log_alpha=log(1e-15),log_beta=log(1e-15),
    log_measure=log(mm))  
}
if(gmodel==2){
  if(vb_log_k_parms$fix_beta==T) map<-append(map,list(log_beta=factor(NA)))
  if(measurer_error_start_bounds$use_measurer==F){
    map<-append(map,list(log_measure=factor(NA)))
    mm<-1e-15
  }
  if(measurer_error_start_bounds$use_measurer==T){
    mm<-measurer_error_start_bounds$sigma_measure
  }
  
  parameters<-list(
             log_K1=log(K_start_bounds$K1),
             log_K2=log(K_start_bounds$K2),
             log_mu_Linf=log(mu_Linf_start_bounds$mu),
             log_mean_age=log(A_start_bounds$mean_age),
             log_sigma_mu_Linf=log(mu_Linf_start_bounds$sigma_mu_Linf),
             log_sigma_age=log(A_start_bounds$sigma_age),
             log_sigma_resid=log(resid_error_start_bounds$sigma_resid),
             A=rep(A_start_bounds$mean_age,length(l1)),#age random effects
             log_alpha=log(vb_log_k_parms$alpha),log_beta=log(vb_log_k_parms$beta),
             log_measure=log(mm))  
}

#Set parameter boundaries
lower_parms<-log(c(K_start_bounds$lower_K1,K_start_bounds$lower_K2,mu_Linf_start_bounds$lower_mu_Linf,A_start_bounds$lower_mean_age,mu_Linf_start_bounds$lower_sigma_mu_Linf,
               A_start_bounds$lower_sigma_age,resid_error_start_bounds$lower_sigma_resid,vb_log_k_parms$lower_alpha,vb_log_k_parms$lower_beta,
               measurer_error_start_bounds$lower_sigma_measure))

lower_parms[lower_parms==-Inf]<-log(1e-15)
upper_parms<-log(c(K_start_bounds$upper_K1,K_start_bounds$upper_K2,mu_Linf_start_bounds$upper_mu_Linf,A_start_bounds$upper_mean_age,mu_Linf_start_bounds$upper_sigma_mu_Linf,
               A_start_bounds$upper_sigma_age,resid_error_start_bounds$upper_sigma_resid,vb_log_k_parms$upper_alpha,vb_log_k_parms$upper_beta,
               measurer_error_start_bounds$upper_sigma_measure))

if(length(map)==0) obj <- MakeADFun(data=c(model="grow_LEP",data),  parameters=parameters,random="A", DLL="fishmethods_TMBExports",inner.control=tmb.control,silent=T)
if(length(map)>0) obj <- MakeADFun(data=c(model="grow_LEP",data),parameters=parameters,random="A",map=map, DLL="fishmethods_TMBExports",inner.control=tmb.control,silent=T)

if(use_parameter_boundaries==T) opt <- nlminb(obj$par, obj$fn, obj$gr,control=nlminb.control,lower=lower_parms,upper=upper_parms)
if(use_parameter_boundaries==F) opt <- nlminb(obj$par, obj$fn, obj$gr,control=nlminb.control)


## Get parameter uncertainties and convergence diagnostics
output <- sdreport(obj,getReportCovariance=T)
outparms<-summary(output)

#Growth Parameters non-log-scale
index1<-which(rownames(outparms)=="mu_Linf")
indexlast<-length(rownames(outparms))
growthparms<-outparms[index1:indexlast,]

#Substitute mean_age with log_mean_age to duplicate Lasett et al. table
index1<-which(rownames(outparms)=="log_mean_age")
lmuage<-outparms[index1,]
index1<-which(rownames(growthparms)=="mean_age")
rownames(growthparms)[index1]<-"mu_log_A"
growthparms[index1,]<-lmuage
index1<-which(rownames(growthparms)=="sigma_age")
rownames(growthparms)[index1]<-"sigma_log_A"

#Insert beta value if fixing beta
  if(gmodel==2 & vb_log_k_parms$fix_beta==T){
    index1<-which(rownames(growthparms)=="alpha")
    d1<-growthparms[1:index1,]
    d2<-growthparms[c((index1+1):nrow(growthparms)),]
    d1<-rbind(d1,cbind(beta=vb_log_k_parms$beta,sd=NA))
    rownames(d1)[nrow(d1)]<-"beta"
    growthparms<-rbind(d1,d2)
  }

# Age Estimates
index1<-which(rownames(outparms)=="A")
As<-outparms[index1,]

#Predicted and residuals
orig_predicted<-data.frame(cbind(pred_l1=obj$report()$mu1,pred_l2=obj$report()$mu2,pred_mu=obj$report()$pred_mu))
orig_residuals<-data.frame(cbind(resid_l1=obj$report()$lresids1,resid_l2=obj$report()$lresids2))


#Estimation Performance
AIC<-function(opt){ #modified Thorson function from TMBhelper
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = 2*k + 2*negloglike
  return(Return)
}
outs<-matrix("NA",nrow=6,ncol=1)
dimnames(outs)<-list(c("-log-likelihood","convergence","iterations","function_eval","function_gradient","message"))
outs[1,1]<-round(opt$objective,5);outs[2,1]<-opt$convergence;outs[3,1]<-opt$iterations;outs[4,1]<-as.numeric(opt$evaluations[1]);
outs[5,1]<-as.numeric(opt$evaluations[2]);outs[6,1]<-opt$message
outs<-data.frame(outs)
names(outs)<-"value"


# Code for generating Atilde provided by Paige Eveson, CSIRO, Hobart, AU
# A as an unknown fixed effect and the other parameters as known.
# This estimator is a close approximation to the conditionally unbiased 
# estimator atilde of A, which can be justified as the correct estimator of A
# for graphical purposes. However, aml is much faster to calculate than Atilde.  

# specify growth functions 
growthvb.f<- function(a,param.g)
{   # this is the VB growth curve
  # generate the parameters
  
  # a is nf x n matrix
  tf<- a<0
  a[tf]<- 0
  
  k<- param.g[1]
  g<- 1-exp(-k*a)
  return(g)
}

growthvblogk.f<- function(a,param.g)
{   # this is the VB log k growth curve
  # generate the parameters
  
  # a is nf x n matrix
  tf<- a<0
  a[tf]<- 0
  k1<-    param.g[1]
  k2<-    param.g[2]
  alpha<- param.g[3]
  beta<-  param.g[4]
  
  theta<- -(k2-k1)/beta
  d1<- 1+exp(-beta*(a-alpha))
  d2<- 1+exp(beta*alpha)
  g<-  (d1/d2)^theta
  g<-  1-exp(-k2*a)*g
  return(g)
}

#
# loglla.f calculates the log-likelihood of A given the other parameters
#
loglla.f<- function(a,l1,l2,dt,mu.Linf,sig.Linf,param.g,sigma)
{        
  a1<- a
  tf<- a<0
  a1[tf]<- 0
  a2<- a1 + dt
  f1<- growth.f(a1, param.g)
  f2<- growth.f(a2, param.g)
  mu1a<- mu.Linf*f1
  mu2a<- mu.Linf*f2
  q1<- l1-mu1a
  q2<- l2-mu2a
  detv<- (sigma*f2)^2+(sigma*f1)^2
  detv<- (sigma*sigma)^2+(sig.Linf^2)*detv
  logh<- (q1*f2-q2*f1)^2
  logh<- logh*(sig.Linf^2)
  logh<- (sigma*q2)^2+(sigma*q1)^2+logh
  logh<- logh/detv
  logh<- -0.5*log(detv)-0.5*logh
  
  # logh is the log likelihood
  # it is clear that it is negative component-wise,
  # so we return sqrt(-logh) which can be treated
  # as a residual sum-of-squares to be minimised by nlminb
  return(sqrt(-logh))
}


calc.atilde.f <- function(l1, l2, dt, mu.Linf, sig.Linf, param.g, sigma)
{
  # this function calculates the approximately conditionally unbiased release 
  # age estimates described in Laslett et al (2004) 
  
  # find an approximate maximum of the log-likelihood for each fish
  # this is given by amax
  #
  aset<- seq(1.0,5.0,0.01)
  amat<- rep(1,length(l1))%*%t(aset)
  llmat<- loglla.f(amat,l1,l2,dt,mu.Linf,sig.Linf,param.g,sigma)
  llmin<- apply(llmat,1,min)
  eps<- .Machine$double.eps
  tf<- (llmat<llmin+eps) | (llmat==llmin)
  
  #
  # row and column numbers of the maxima
  #
  nc<- col(tf)[tf]
  nr<- row(tf)[tf]
  ord<- order(nr,nc)
  nr<- nr[ord]
  nc<- nc[ord]
  
  #
  # eliminate duplicate maxima within rows (if any)
  #
  tf<- !duplicated(nr)
  nr<- nr[tf]
  nc<- nc[tf]
  amax<- amat[cbind(nr,nc)]
  amax<- as.vector(amax)
  
  #
  # calculate the maximum likelihood estimators of A using nlminb
  #
  m<- length(dt)
  aml<- numeric(0)
  mess<- character(0)
 
  for(i in (1:m) )
  {
    # using nlminb for optimisation
    nl.fit<- nlminb(start = amax[i], objective = loglla.f, 
                    dt=dt[i], l1=l1[i], l2=l2[i], mu.Linf=mu.Linf, sig.Linf=sig.Linf,
                    param.g=param.g, sigma=sigma, lower=0.5, upper=6)
    
    aml<- c(aml,nl.fit$par)
    mess<- c(mess,nl.fit$mess)
    lb<- paste(i," out of ",m," completed ")
    lb<- paste(lb,"     ",nl.fit$mess)
    print(lb)
  }
  print("Atilde Computations Completed.")
  return(aml)
}

if(gmodel == 1){
  growth.f <- growthvb.f
  param.g<-growthparms[rownames(growthparms)=="K1",1]
  mu.Linf <-  growthparms[rownames(growthparms)=="mu_Linf",1]
  sig.Linf <-growthparms[rownames(growthparms)=="sigma_mu_Linf",1]
  sigma<-growthparms[rownames(growthparms)=="sigma_resid",1]
  atilde <- calc.atilde.f(l1, l2, dt, mu.Linf, sig.Linf, param.g, sigma)
  atilde_lpred1=mu.Linf*(1-exp(-param.g*atilde))
  atilde_lresid1=c(l1-atilde_lpred1)
  atilde_lpred2<-mu.Linf*(1-exp(-param.g*c(atilde+dt)))
  atilde_lresid2=c(l2-atilde_lpred2)
  atilde_results<-data.frame(atilde=atilde,atilde_predicted_l1=atilde_lpred1,atilde_predicted_l2=atilde_lpred2,atilde_residuals_l1=atilde_lresid1,atilde_residuals_l2=atilde_lresid2)
 
  if(graphs){ 
    par(mfrow=c(2,2))
    plot(l1~atilde,ylab="l1",xlab=bquote({tilde(A)}['f']))
    pred1<-data.frame(cbind(lpred1=atilde_lpred1,atilde=atilde))
    pred1<-pred1[order(pred1$atilde),]
    lines(pred1$lpred1~pred1$atilde)
    plot(atilde_lresid1~atilde_lpred1,ylab="Residuals",xlab="Fitted Length")
    abline(h=0,col="red")

    plot(l2~c(atilde+dt),ylab="l2",xlab=bquote({tilde(A)}['f']+{dt}))
    pred2<-data.frame(cbind(lpred2=atilde_lpred2,atilde=c(atilde+dt)))
    pred2<-pred2[order(pred2$atilde),]
    lines(pred2$lpred2~pred2$atilde) 
    plot(atilde_lresid2~atilde_lpred2,ylab="Residuals",xlab="Fitted Length")
    abline(h=0,col="red")
  }
}
if(gmodel == 2){
  growth.f <- growthvblogk.f
  param.g<-c(growthparms[rownames(growthparms)=="K1",1],growthparms[rownames(growthparms)=="K2",1],growthparms[rownames(growthparms)=="alpha",1],
             growthparms[rownames(growthparms)=="beta",1])
  mu.Linf <-  growthparms[rownames(growthparms)=="mu_Linf",1]
  sig.Linf <-growthparms[rownames(growthparms)=="sigma_mu_Linf",1]
  sigma<-growthparms[rownames(growthparms)=="sigma_resid",1] 
  atilde <- calc.atilde.f(l1, l2, dt, mu.Linf, sig.Linf, param.g, sigma)
  eq11=(1+exp(-param.g[4]*(atilde-param.g[3])))/(1+exp(param.g[4]*param.g[3]));
  eq12=(-(param.g[2]-param.g[1]))/param.g[4];
  eq21=((1+exp(-param.g[4]*(atilde+dt-param.g[3])))/(1+exp(param.g[4]*param.g[3])));
  
  atilde_lpred1=mu.Linf*(1-exp(-param.g[2]*(atilde))*(eq11^eq12))
  atilde_lpred2=mu.Linf*(1-exp(-param.g[2]*(atilde+dt))*(eq21^eq12))
  atilde_lresid1=c(l1-atilde_lpred1)
  atilde_lresid2=c(l2-atilde_lpred2)
  atilde_results<-data.frame(atilde=atilde,atilde_predicted_l1=atilde_lpred1,atilde_predicted_l2=atilde_lpred2,atilde_residuals_l1=atilde_lresid1,atilde_residuals_l2=atilde_lresid2)
  
  if(graphs){
  par(mfrow=c(2,2))
  plot(l1~atilde,ylab="l1",xlab=bquote({tilde(A)}['f']))
  pred1<-data.frame(cbind(lpred1=atilde_lpred1,atilde=atilde))
  pred1<-pred1[order(pred1$atilde),]
  lines(pred1$lpred1~pred1$atilde)
  plot(atilde_lresid1~atilde_lpred1,ylab="Residuals",xlab="Fitted Length")
  abline(h=0,col="red")
  
  plot(l2~c(atilde+dt),ylab="l2",xlab=bquote({tilde(A)}['f']+{dt}))
  pred2<-data.frame(cbind(lpred2=atilde_lpred2,atilde=c(atilde+dt)))
  pred2<-pred2[order(pred2$atilde),]
  lines(pred2$lpred2~pred2$atilde)
  plot(atilde_lresid2~atilde_lpred1,ylab="Residuals",xlab="Fitted Length")
  abline(h=0,col="red")
  }
}
outpt_list<-list(parameter_estimates=growthparms,AIC=AIC(opt),A=As,orig_predicted=orig_predicted,orig_residuals=orig_residuals,atilde_results=atilde_results,converge_stats=outs)

return(outpt_list)
}#function
