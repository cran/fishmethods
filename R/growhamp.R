growhamp<-function(L1 = NULL, L2 = NULL, TAL = NULL, models=c(1,2,3,4,5,6,7),
          method=c("Nelder-Mead","Nelder-Mead","Nelder-Mead","Nelder-Mead","Nelder-Mead",
                   "Nelder-Mead","Nelder-Mead"),
          varcov=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
          Linf=list(startLinf=NULL,lowerLinf=NULL,upperLinf=NULL),       
 	    K=list(startK=NULL,lowerK=NULL,upperK=NULL),
	    sigma2_error=list(startsigma2=NULL,lowersigma2=NULL,uppersigma2=NULL),
	    sigma2_Linf=list(startsigma2=NULL,lowersigma2=NULL,uppersigma2=NULL),	
	    sigma2_K=list(startsigma2=NULL,lowersigma2=NULL,uppersigma2=NULL),
          mu_measure=0, sigma2_measure=0,
          control=list(maxit=1000)){
# Check models and methods
if(!all(method %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN"))) stop("method doesn't exist") 
if(any(models<1)||any(models>7)) stop("model doesn't exist") 
if(any(duplicated(models))) stop("duplicates exist in 'models'")
  
ml<-c(length(models),length(method))
if(!isTRUE(all.equal(ml,rep(ml[1],length(ml))))) warning("Lengths of 'method' and 'models' do not match")
if(length(method)<length(models)) method<-as.character(c(method,rep(method[length(method)],c(length(models)-length(method)))))
if(length(method)>length(models)) method<-as.character(method[1:length(models)])
ml<-c(length(models),length(varcov))
if(!isTRUE(all.equal(ml,rep(ml[1],length(ml))))) warning("Lengths of 'varcov' and 'models' do not match")
if(length(varcov)<length(models)) varcov<-c(varcov,rep(varcov[length(varcov)],c(length(models)-length(varcov))))
if(length(varcov)>length(models)) varcov<-varcov[1:length(models)]
mm<-data.frame(model=models,method=method)
mm$method<-as.character(mm$method)
mm<-mm[order(mm$model),]

# Check for equal lengths
cl<-c(length(L1),length(L2),length(TAL))
if(!isTRUE(all.equal(cl,rep(cl[1],length(cl))))) stop("Lengths of observation vectors do not match")
grdat<-data.frame(L1=L1,L2=L2,ti=TAL)
#Check for missing data
if(nrow(grdat)!=nrow(na.omit(grdat))){
dr<-nrow(grdat)-nrow(na.omit(grdat))
warning(paste(dr," row(s) deleted due to missingness",sep=""))
grdat<-na.omit(grdat)
}
#Warning about small sample sizes
if(nrow(grdat)<10) warning("Less than 10 rows of data are available")
grdat$dl<-grdat[,2]-grdat[,1]
n=nrow(grdat)
results <- data.frame(Model = models, Linf= NA,
       K = NA, s2Linf= NA,s2K=NA, 
       s2error= NA,boundary=NA,LL= NA,AIC=NA,method="NA",stringsAsFactors = FALSE)
names(results)[8]<-"-Log Likelihood"
residuals<-matrix(0,nrow=n,ncol=length(models))
varcov1<-NULL
cnt=0
#Models
if(any(mm$model==1)){
cnt=cnt+1
warn1=0
lab1="Faber"
parms<-c(Linf$startLinf,K$startK,sigma2_error$startsigma2)
svb<-function(parms){
pred<-(parms[1]-grdat$L1)*(1-exp(-parms[2]*grdat$ti))
p1<-log(2*pi*parms[3])
(n/2)*log(2*pi*parms[3])+sum((grdat$dl-pred)^2)/(2*parms[3])
}
index<-which(mm$model==1)
if(mm$method[index] %in% c("Nelder-Mead","BFGS","CG","SANN")) model1<-optim(parms,svb,method=mm$method[index],hessian=varcov[cnt],control=control)
if(mm$method[index] %in% c("L-BFGS-B")){
   model1<-optim(parms,svb,method=mm$method[index],hessian=varcov[cnt],lower=c(Linf$lowerLinf,
   K$lowerK,sigma2_error$lowersigma2),upper=c(Linf$upperLinf,
   K$upperK,sigma2_error$uppersigma2),control=control)
   bound<-c(model1$par[1]==Linf$lowerLinf,model1$par[1]==Linf$upperLinf,
            model1$par[2]==K$lowerK,model1$par[2]==K$upperK,
            model1$par[3]==sigma2_error$lowersigma2,model1$par[3]==sigma2_error$uppersigma2)  
    if(any(bound==TRUE)) warn1<-1 else warn1<-0
}
results[cnt,1]<-lab1
results[cnt,2]<-model1$par[1];results[cnt,3]<-model1$par[2];
results[cnt,6]<-model1$par[3];results[cnt,7]<-warn1;results[cnt,8]<-model1$value;
results[cnt,9]<-2*model1$value+2*length(model1$par)
results[cnt,10]<-mm$method[index]
residuals[,cnt]<-grdat$dl-(model1$par[1]-grdat$L1)*(1-exp(-model1$par[2]*grdat$ti))
if(varcov[cnt]==TRUE) varcov1[cnt]<-list(solve(model1$hessian))
}#model 1

if(any(mm$model==2)){
cnt=cnt+1
warn2<-0
lab2="Kirkwood and Somers"
parms<-c(Linf$startLinf,K$startK,sigma2_Linf$startsigma2)
ksm<-function(parms){
Edl<-(parms[1]-grdat$L1)*(1-exp(-parms[2]*grdat$ti))
vardl<-parms[3]*(1-exp(-parms[2]*grdat$ti))^2
vardl[vardl<=0]<-1e-8
p1<-log(2*pi*vardl)/2
p2<-((grdat$dl-Edl)^2)/(2*vardl)
sum(p1+p2)
}

index<-which(mm$model==2)
if(mm$method[index] %in% c("Nelder-Mead","BFGS","CG","SANN")) model2<-optim(parms,ksm,method=mm$method[index],hessian=varcov[cnt],control=control)
if(mm$method[index] %in% c("L-BFGS-B")){
    model2<-optim(parms,ksm,method=mm$method[index],hessian=varcov[cnt],lower=c(Linf$lowerLinf,
    K$lowerK,sigma2_Linf$lowersigma2),upper=c(Linf$upperLinf,
   K$upperK,sigma2_Linf$uppersigma2),control=control)
    bound<-c(model2$par[1]==Linf$lowerLinf,model2$par[1]==Linf$upperLinf,
             model2$par[2]==K$lowerK,model2$par[2]==K$upperK,
             model2$par[3]==sigma2_Linf$lowersigma2,model2$par[3]==sigma2_Linf$uppersigma2)
    if(any(bound==TRUE)) warn2<-1 else warn2<-0
}
results[cnt,1]<-lab2
results[cnt,2]<-model2$par[1];results[cnt,3]<-model2$par[2];
results[cnt,4]<-model2$par[3];
results[cnt,7]<-warn2;results[cnt,8]<-model2$value;
results[cnt,9]<-2*model2$value+2*length(model2$par)
results[cnt,10]<-mm$method[index]
residuals[,cnt]<-grdat$dl-(model2$par[1]-grdat$L1)*(1-exp(-model2$par[2]*grdat$ti))
if(varcov[cnt]==TRUE) varcov1[cnt]<-list(solve(model2$hessian))
}#model 2

if(any(mm$model==3)){
cnt=cnt+1
warn3<-0
lab3="Kirkwood and Somers with ME"
parms<-c(Linf$startLinf,K$startK,sigma2_Linf$startsigma2,sigma2_error$startsigma2)
ksme<-function(parms){
Edl<-(parms[1]-grdat$L1)*(1-exp(-parms[2]*grdat$ti))
vardl<-parms[3]*((1-exp(-parms[2]*grdat$ti))^2)+parms[4]
vardl[vardl<=0]<-1e-8
p1<-log(2*pi*vardl)/2
p2<-((grdat$dl-Edl)^2)/(2*vardl)
sum(p1+p2)
}
index<-which(mm$model==3)
if(mm$method[index] %in% c("Nelder-Mead","BFGS","CG","SANN")) model3<-optim(parms,ksme,hessian=varcov[cnt],method=mm$method[index],control=control)
if(mm$method[index] %in% c("L-BFGS-B")){
   model3<-optim(parms,ksme,method=mm$method[index],hessian=varcov[cnt],lower=c(Linf$lowerLinf,
   K$lowerK,sigma2_Linf$lowersigma2,sigma2_error$lowersigma2),
   upper=c(Linf$upperLinf,K$upperK,sigma2_Linf$uppersigma2,sigma2_error$uppersigma2),control=control)
   bound<-c(model3$par[1]==Linf$lowerLinf,model3$par[1]==Linf$upperLinf,
            model3$par[2]==K$lowerK,model3$par[2]==K$upperK,
            model3$par[3]==sigma2_Linf$lowersigma2,model3$par[3]==sigma2_Linf$uppersigma2,
            model3$par[4]==sigma2_error$lowersigma2,model3$par[4]==sigma2_error$uppersigma2)
   if(any(bound==TRUE)) warn3<-1 else warn3<-0
}
results[cnt,1]<-lab3
results[cnt,2]<-model3$par[1];results[cnt,3]<-model3$par[2];
results[cnt,4]<-model3$par[3];results[cnt,6]<-model3$par[4]
results[cnt,7]<-warn3;results[cnt,8]<-model3$value;
results[cnt,9]<-2*model3$value+2*length(model3$par)
results[cnt,10]<-mm$method[index]
residuals[,cnt]<-grdat$dl-(model3$par[1]-grdat$L1)*(1-exp(-model3$par[2]*grdat$ti))
if(varcov[cnt]==TRUE) varcov1[cnt]<-list(solve(model3$hessian))
}#model 3


if(any(mm$model==4)){
cnt=cnt+1
warn4<-0
lab4="Kirkwood and Somers with ME & RLE"
parms<-c(Linf$startLinf,K$startK,sigma2_Linf$startsigma2,sigma2_error$startsigma2)
ksme<-function(parms){
EL2<-(parms[1]-grdat$L1)*(1-exp(-parms[2]*grdat$ti))+grdat$L1+exp(-parms[2]*grdat$ti)*mu_measure
varL2<-parms[3]*(1-exp(-parms[2]*grdat$ti))^2+parms[4]+sigma2_measure*exp(-2*parms[2]*grdat$ti)
varL2[varL2<=0]<-1e-8
p1<-log(2*pi*varL2)/2
p2<-((grdat$L2-EL2)^2)/(2*varL2)
sum(p1+p2)
}

index<-which(mm$model==4)
if(mm$method[index] %in% c("Nelder-Mead","BFGS","CG","SANN")) model4<-optim(parms,ksme,hessian=varcov[cnt],method=mm$method[index],control=control)
if(mm$method[index] %in% c("L-BFGS-B")){
  model4<-optim(parms,ksme,method=mm$method[index],hessian=varcov[cnt],lower=c(Linf$lowerLinf,
   K$lowerK,sigma2_Linf$lowersigma2,sigma2_error$lowersigma2),
   upper=c(Linf$upperLinf,K$upperK,sigma2_Linf$uppersigma2,sigma2_error$uppersigma2),control=control)
  bound<-c(model4$par[1]==Linf$lowerLinf,model4$par[1]==Linf$upperLinf,
           model4$par[2]==K$lowerK,model4$par[2]==K$upperK,
           model4$par[3]==sigma2_Linf$lowersigma2,model4$par[3]==sigma2_Linf$uppersigma2,
           model4$par[4]==sigma2_error$lowersigma2,model4$par[4]==sigma2_error$uppersigma2
           )
  
  
   if(any(bound==TRUE)) warn4<-1 else warn4<-0
}
results[cnt,1]<-lab4
results[cnt,2]<-model4$par[1];results[cnt,3]<-model4$par[2];
results[cnt,4]<-model4$par[3];results[cnt,6]<-model4$par[4]
results[cnt,7]<-warn4;results[cnt,8]<-model4$value;
results[cnt,9]<-2*model4$value+2*length(model4$par)
results[cnt,10]<-mm$method[index]
residuals[,cnt]<-grdat$L2-((model4$par[1]-grdat$L1)*(1-exp(-model4$par[2]*grdat$ti))+grdat$L1+exp(-model4$par[2]*grdat$ti)*mu_measure+model4$par[4])
if(varcov[cnt]==TRUE) varcov1[cnt]<-list(solve(model4$hessian))
}#model 4


if(any(mm$model==5)){
cnt=cnt+1
warn5<-0
lab5="Sainsbury"
parms<-c(Linf$startLinf,K$startK,sigma2_Linf$startsigma2,sigma2_K$startsigma2)
sain<-function(parms){
ep<-(-(parms[2]^2)/parms[4])
ip<-1+(parms[4]*grdat$ti/parms[2])
ip2<-1+(2*parms[4]*grdat$ti/parms[2])
Edl<-(parms[1]-grdat$L1)*(1-(ip^ep))
C1<-1-2*(ip^ep)+ip2^ep

C2<-(ip2^ep)-ip^(2*ep)
vardl<-C1*parms[3]+C2*(parms[1]-grdat$L1)^2
vardl[vardl<=0]<-1e-8
p1<-log(2*pi*vardl)/2
p2<-((grdat$dl-Edl)^2)/(2*vardl)
sum(p1+p2)
}

index<-which(mm$model==5)
if(mm$method[index] %in% c("Nelder-Mead","BFGS","CG","SANN")) model5<-optim(parms,sain,hessian=varcov[cnt],method=mm$method[index],control=control)
if(mm$method[index] %in% c("L-BFGS-B")){
   model5<-optim(parms,sain,method=mm$method[index],hessian=varcov[cnt],lower=c(Linf$lowerLinf,
   K$lowerK,sigma2_Linf$lowersigma2,sigma2_K$lowersigma2),
   upper=c(Linf$upperLinf,K$upperK,sigma2_Linf$uppersigma2,sigma2_K$uppersigma2),control=control)
  bound<-c(model5$par[1]==Linf$lowerLinf,model5$par[1]==Linf$upperLinf,
                     model5$par[2]==K$lowerK,model5$par[2]==K$upperK,
                     model5$par[3]==sigma2_Linf$lowersigma2,model5$par[3]==sigma2_Linf$uppersigma2,
                     model5$par[4]==sigma2_K$lowersigma2,model5$par[4]==sigma2_K$uppersigma2)
  
   if(any(bound==TRUE)) warn5<-1 else warn5<-0
}
results[cnt,1]<-lab5
results[cnt,2]<-model5$par[1];results[cnt,3]<-model5$par[2];
results[cnt,4]<-model5$par[3];results[cnt,5]<-model5$par[4]
results[cnt,7]<-warn5;results[cnt,8]<-model5$value;
results[cnt,9]<-2*model5$value+2*length(model5$par)
results[cnt,10]<-mm$method[index]
ip<-1+(model5$par[4]*grdat$ti/model5$par[2])
ep<-(-(model5$par[2]^2)/model5$par[4])
Edl<-(model5$par[1]-grdat$L1)*(1-(ip^ep))
residuals[,cnt]<-grdat$dl-Edl
if(varcov[cnt]==TRUE) varcov1[cnt]<-list(solve(model5$hessian))
}#model 5

if(any(mm$model==6)){
cnt=cnt+1
warn6=0
lab6="Sainsbury with ME"
parms<-c(Linf$startLinf,K$startK,sigma2_Linf$startsigma2,sigma2_K$startsigma2,sigma2_error$startsigma2)
sain<-function(parms){
ep<--(parms[2]^2)/parms[4]
ip<-1+(parms[4]*grdat$ti/parms[2])
ip2<-1+(2*parms[4]*grdat$ti/parms[2])
Edl<-(parms[1]-grdat$L1)*(1-(ip^ep))
C1=1-2*(ip^ep)+ip2^ep
C2=(ip2^ep)-ip^(2*ep)
vardl<-C1*parms[3]+(C2*(parms[1]-grdat$L1)^2)+parms[5]
vardl[vardl<=0]<-1e-8
p1<-log(2*pi*vardl)/2
p2<-((grdat$dl-Edl)^2)/(2*vardl)
sum(p1+p2)
}

index<-which(mm$model==6)
if(mm$method[index] %in% c("Nelder-Mead","BFGS","CG","SANN")) model6<-optim(parms,sain,hessian=varcov[cnt],method=mm$method[index],control=control)
if(mm$method[index] %in% c("L-BFGS-B")){
   model6<-optim(parms,sain,method=mm$method[index],hessian=varcov[cnt],lower=c(Linf$lowerLinf,
   K$lowerK,sigma2_Linf$lowersigma2,sigma2_K$lowersigma2,sigma2_error$lowersigma2),
   upper=c(Linf$upperLinf,K$upperK,sigma2_Linf$uppersigma2,sigma2_K$uppersigma2,sigma2_error$uppersigma2)
   ,control=control)
            bound<-c(model6$par[1]==Linf$lowerLinf,model6$par[1]==Linf$upperLinf,
                     model6$par[2]==K$lowerK,model6$par[2]==K$upperK,
                     model6$par[3]==sigma2_Linf$lowersigma2,model6$par[3]==sigma2_Linf$uppersigma2,
                     model6$par[4]==sigma2_K$lowersigma2,model6$par[4]==sigma2_K$uppersigma2,
                     model6$par[5]==sigma2_error$lowersigma2,model6$par[5]==sigma2_error$uppersigma2)
            
   if(any(bound==TRUE)) warn6<-1 else warn6<-0
}
results[cnt,1]<-lab6
results[cnt,2]<-model6$par[1];results[cnt,3]<-model6$par[2];
results[cnt,4]<-model6$par[3];results[cnt,5]<-model6$par[4];
results[cnt,6]<-model6$par[5]
results[cnt,7]<-warn6;results[cnt,8]<-model6$value;
results[cnt,9]<-2*model6$value+2*length(model6$par)
results[cnt,10]<-mm$method[index]
ip<-1+(model6$par[4]*grdat$ti/model6$par[2])
ep<-(-(model6$par[2]^2)/model6$par[4])
Edl<-(model6$par[1]-grdat$L1)*(1-(ip^ep))
residuals[,cnt]<-grdat$dl-Edl
if(varcov[cnt]==TRUE) varcov1[cnt]<-list(solve(model6$hessian))
}#model 6

if(any(mm$model==7)){
cnt=cnt+1
warn7<-0
lab7="Sainsbury with ME & RLE"
parms<-c(Linf$startLinf,K$startK,sigma2_Linf$startsigma2,sigma2_K$startsigma2,sigma2_error$startsigma2)
ksme<-function(parms){
EL1<-grdat$L1+mu_measure
ep<--(parms[2]^2)/parms[4]
ip<-1+(parms[4]*grdat$ti/parms[2])
ip2<-1+(2*parms[4]*grdat$ti/parms[2])
EL2<-(parms[1]-EL1)*(1-(ip^ep))+EL1
Ek<-ip2^ep
C1=1-2*(ip^ep)+(ip2^ep)
C2=(ip2^ep)-(ip^(2*ep))
varL2<-C1*parms[3]+(C2*(parms[1]-EL1)^2)+parms[5]+sigma2_measure*Ek
varL2[varL2<=0]<-1e-8
p1<-log(2*pi*varL2)/2
p2<-((grdat$L2-EL2)^2)/(2*varL2)
sum(p1+p2)
}

index<-which(mm$model==7)
if(mm$method[index] %in% c("Nelder-Mead","BFGS","CG","SANN")) model7<-optim(parms,ksme,hessian=varcov[cnt],method=mm$method[index],control=control)
if(mm$method[index] %in% c("L-BFGS-B")){
   model7<-optim(parms,ksme,method=mm$method[index],hessian=varcov[cnt],lower=c(Linf$lowerLinf,
   K$lowerK,sigma2_Linf$lowersigma2,sigma2_K$lowersigma2,sigma2_error$lowersigma2),
   upper=c(Linf$upperLinf,K$upperK,sigma2_Linf$uppersigma2,sigma2_K$uppersigma2,sigma2_error$uppersigma2),control=control)
bound<-c(model7$par[1]==Linf$lowerLinf,model7$par[1]==Linf$upperLinf,
         model7$par[2]==K$lowerK,model7$par[2]==K$upperK,
         model7$par[3]==sigma2_Linf$lowersigma2,model7$par[3]==sigma2_Linf$uppersigma2,
         model7$par[4]==sigma2_K$lowersigma2,model7$par[4]==sigma2_K$uppersigma2,
         model7$par[5]==sigma2_error$lowersigma2,model7$par[5]==sigma2_error$uppersigma2)
   if(any(bound==TRUE)) warn7<-1 else warn7<-0
}
   results[cnt,1]<-lab7
   results[cnt,2]<-model7$par[1];results[cnt,3]<-model7$par[2];
   results[cnt,4]<-model7$par[3];results[cnt,5]<-model7$par[4];
   results[cnt,6]<-model7$par[5]
   results[cnt,7]<-warn7;results[cnt,8]<-model7$value;
   results[cnt,9]<-2*model7$value+2*length(model7$par)
   results[cnt,10]<-mm$method[index]
ip<-1+(model7$par[4]*grdat$ti/model7$par[2])
ep<-(-(model7$par[2]^2)/model7$par[4])
EL2<-(model7$par[1]-(grdat$L1+mu_measure))*(1-(ip^ep))+(grdat$L1+mu_measure)
residuals[,cnt]<-grdat$L2-EL2
if(varcov[cnt]==TRUE) varcov1[cnt]<-list(solve(model7$hessian))
}#model 7
labels<-data.frame(models=models,varcov=varcov)

if(!is.null(varcov1)) names(varcov1)<-c(paste("model",labels[labels$varcov=="TRUE",1]))
residuals<-as.data.frame(residuals)
names(residuals)<-paste("model",models)
outpt<-list(results=results,varcov=varcov1,residuals=residuals)
return(outpt)
}
