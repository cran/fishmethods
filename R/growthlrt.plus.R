
growthlrt.plus<-function(model,data=sys.frame(sys.parent()),
                         params=NULL,
                         start=NULL,
                         within_grp_var=~1,
                         cfh=NULL,
                         nlminb.control=list(iter.max=10000,rel.tol=1e-10),
                         optim.control=list(maxit=10000,reltol=1e-10)){
  formula<-as.formula(model)
  if(!inherits(formula,"formula")) stop("'model' is not a formula")
  formula<-as.formula(formula)
  
  if(length(formula)!=3) stop("object model must be of the form \"resp~pred\"")
  if(is.null(start)) stop("no start values for model parameters")
  if(!is.list(start)) stop("start values must be a list")
  
  if(is.null(params)) stop("no params formulae") 
  if(!is.list(params)) stop("params must be a list")
  if(!is.data.frame(data)) stop("'data' is not a data.frame (or list)")
  
  ###----------------------- check params argument ------------------------------##
  pnames <- character(length(params))
  for (i in seq_along(params)) {
    this <- eval(params[[i]])
    if (!inherits(this, "formula")) 
      stop("'params' must be a formula or list of formulae")
    if (length(this) != 3) 
      stop("formulae in 'params' must be of the form \"parameter ~ 1\"")
    if (!is.name(this[[2]])) 
      stop("formulae in 'params' must be of the form \"parameter ~ 1\"")
    pnames[i] <- as.character(this[[2]])
  }
  names(params) <- pnames
  if(!all(names(start)==pnames)) stop("parameter names in 'start' and 'params' do not match")
  
  grp_name<-NULL
  ###-----------------check within_grp_var --------------------------------------------------###
  wg_formula<-as.formula(within_grp_var)
  if(!inherits(wg_formula,"formula")) stop("'within_grp_var' is not a formula")
  if(length(wg_formula)!=2) stop("object within_grp_var must be of the form \"~1 or ~grp\"")
  if(length(as.character(wg_formula[[2]]))!=1) stop("object within_grp_var must be of the form \"~1 or ~grp\"")
  wgnames<-as.character(wg_formula[[2]])
  if(wgnames=="1"){nsigmas<-1}
  if(wgnames!="1"){
    grp_name<-wgnames
    #check grp is in data file
    datacolnames<-names(data)
    cumlabels<-NULL
    if(!all(wgnames %in% datacolnames)) cumlabels<-c(wgnames)
    if(!is.null(cumlabels)) stop(paste("within_grp_var: '",paste(cumlabels,collapse=", "),"' not in 'data'",sep="")) 
    nsigmas<-eval(parse(text=paste("nlevels(as.factor(data$",wgnames,"))",sep="")))
    if(nsigmas<2) stop(paste("within_grp_var: only one level present in '",wgnames,"'."," Change formula to ~1",sep=""))
  }
  
  
  ###---------------check cfh ---------------------------------------------------------------###
  if(is.null(cfh)){
    nthetas<-0
    cfh_fixed<-NULL
    cfh_names<-NULL
    cfh_value<-NULL
  }
 
  if(!is.null(cfh)){
    #are required argument present?
    if(!any(names(cfh)=="form")) stop("'cfh form' is not present")
    if(!any(names(cfh)=="value")) stop("'cfh value' is not present")
    if(!any(names(cfh)=="fixed")) stop("'cfh fixed' is not present")
    if(!inherits(cfh$form,"formula")) stop("'cfh formula' is not a formula")
    if(!is.numeric(cfh$value)) stop("'cfh value' is not numeric")
    if(!is.null(cfh$fixed)){
      if(!is.numeric(cfh$fixed)) stop("'cfh fixed' is not numeric or NULL")
      
      }

    if(length(cfh$form)!=2) stop("object cfh formula must be of the form \"~1 or ~grp\"")
    if(length(as.character(cfh$form[[2]]))!=1) stop("object cfh formula must be of the form \"~1 or ~grp\"")
    cfh_names<-as.character(cfh$form[[2]])
    if(cfh_names=="1"){
      nthetas<-1
      cfh_value<-cfh$value
      if(length(cfh$value)>nthetas) cfh_value<-cfh_value[1:nthetas]
      if(length(cfh$value)<nthetas) cfh_value<-rep(cfh_value,nthetas)
      cfh_fixed<-cfh$fixed
      if(!is.null(cfh_fixed)){
        if(length(cfh_fixed)>nthetas) cfh_fixed<-cfh_fixed[1:nthetas]
        if(length(cfh_fixed)<nthetas) cfh_fixed<-rep(cfh_fixed,nthetas)
      }
    }
    if(cfh_names!="1"){
      if(is.null(grp_name)) grp_name<-cfh_names
      #check grp is in data file
      datacolnames<-names(data)
      cumlabels<-NULL
      if(!all(cfh_names %in% datacolnames)) cumlabels<-c(cfh_names)
      if(!is.null(cumlabels)) stop(paste("cfh formula: '",paste(cumlabels,collapse=", "),"' not in 'data'",sep="")) 
      nthetas<-eval(parse(text=paste("nlevels(as.factor(data$",cfh_names,"))",sep="")))
      if(nthetas<2) stop(paste("cfh formula: only one level present in '",cfh_names,"'."," Change formula to ~1",sep=""))
      cfh_value<-cfh$value
      if(length(cfh$value)>nthetas) cfh_value<-cfh_value[1:nthetas]
      if(length(cfh$value)<nthetas) cfh_value<-rep(cfh_value,nthetas)
      cfh_fixed<-cfh$fixed
      if(!is.null(cfh_fixed)){
        if(length(cfh_fixed)>nthetas) cfh_fixed<-cfh_fixed[1:nthetas]
        if(length(cfh_fixed)<nthetas) cfh_fixed<-rep(cfh_fixed,nthetas)
      }
    }
  }
  
  group_name<-NULL
  for(ll in 1:length(params)){
    if(is.numeric(params[[ll]][[3]])){ #if 1 no grouping variable
      if(as.numeric(params[[ll]][[3]])!=1) stop(paste("Parameter ",names(params)[ll], " formula after ~ must be 1 or group name",sep=""))
    }
    if(is.name(params[[ll]][[3]])){  #group variable name
      group_name<-c(group_name,as.character(params[[ll]][[3]]))
      if(is.null(grp_name)) grp_name<-group_name
    }
  }
  if(!is.null(group_name)){
    if(!all(group_name==group_name[1])) stop("Only one group name can be used across parameter formulae")
  }
  
  #Extract names of variable and parameters in rhs equation
  varNames<-all.vars(formula[[3]])
  vars<-varNames[-which(varNames %in% pnames)] #predictors
  
  #If group variable is indicated but there must be a minimum of 2 level
  ngroups<-NULL
  if(!is.null(group_name)){
    group_name<-unique(group_name)
    vars<-c(vars,group_name)
    ngroups<-eval(parse(text=paste("nlevels(as.factor(data$",group_name,"))",sep="")))
    if(ngroups==0) stop(paste("'",group_name,"' not in data",sep=""))
    if(ngroups<2) stop("Only two or more groups are allowed.")
  }
  
  #  Check if the number of start values for each parameter equals the number of parameters that will be estimated
  if(!is.null(ngroups)){
    nparms<-length(params)
    for(j in 1:nparms){
      checkparms<-attr(terms(formula(params[[j]])), which = "term.labels")
      if(length(checkparms)==0){
        ptnames<-names(params)[j]
        tstart<-eval(parse(text=paste("start$",ptnames,sep="")))
        if(length(tstart)>1){
          warning("Number of starting ", ptnames," values > number of ", ptnames, " parameters. Will truncate starting values.")
          newstarts<-tstart[1]
          eval(parse(text=paste("start$",ptnames,"<-newstarts",sep="")))
        }
      }
      if(length(checkparms)>0) {
        ptnames<-names(params)[j]
        tstart<-eval(parse(text=paste("start$",ptnames,sep="")))
        if(length(tstart)<ngroups){
          warning(paste("Number of starting ", ptnames, " values < number of ", ptnames, 
                        " parameters. Will use 1/10th of value as addition.",sep=""))
          newstarts<-c(tstart,rep(tstart[length(tstart)]/10,ngroups-length(tstart)))
          eval(parse(text=paste("start$",ptnames,"<-newstarts",sep="")))
        }
        
        if(length(tstart)>ngroups){
          warning("Number of starting ", ptnames," values > number of Linf parameters. Will truncate starting values.")
          if(!is.null(ngroups)){newstarts<-tstart[1:ngroups]} else 
            eval(parse(text=paste("start$",ptnames,"<-newstarts",sep="")))
        }   
      } 
    }# j loop
  }#!isnull
  
  if(is.null(ngroups)){
    nparms<-length(params)
    for(j in 1:nparms){
      ptnames<-names(params)[j]
      tstart<-eval(parse(text=paste("start$",ptnames,sep="")))
      if(length(tstart)>1){
        warning("Number of starting ", ptnames," values > number of Linf parameters. Will truncate starting values.")
        newstarts<-tstart[1]
        eval(parse(text=paste("start$",ptnames,"<-newstarts",sep="")))
      }
    }# j loop
  }#!isnull
  
  # Is dependent variable name is data set?
  if(length(which(as.character(formula[[2]]) %in% c(names(data))))==0) stop(paste("'",as.character(formula[[2]]),"'",
                                                    " not in data",sep=""))
 
  #Check if predictors from equation are in the dataset
  datacolnames<-names(data)
  cumlabels<-NULL
  for(i in 1:length(vars)){
    if(!all(vars[i] %in% datacolnames)) cumlabels<-c(cumlabels,vars[i])
  }
  if(!is.null(cumlabels)) stop(paste("equation variable(s): ",paste(cumlabels,collapse=", ")," not in 'data'",sep=""))  
  

  #Get list of predictor variables
  if(!is.null(group_name)){
    group_name<-unique(group_name)
    rhs_equation_vars<-vars[-which(vars==group_name)] 
  } else rhs_equation_vars<-vars
  
  #Check all predictors are numeric
  for(j in 1:length(rhs_equation_vars)){
    et<-paste("c(",paste("'",rhs_equation_vars[j],"'",sep="",collapse=","),")",sep="")  
  if(is.factor(data[,eval(parse(text=et))])|is.character(data[,eval(parse(text=et))])) 
       stop(paste("'",rhs_equation_vars[j], "' is not a numeric variable.",sep=""))
  }
  
   #If params, within_grp_var or cfh use a grouping variable ensure it is the same through out
   comb<-NULL   
   if(!is.null(group_name)) comb<-c(comb,group_name)
   if(wgnames!="1") comb<-c(comb,wgnames)
   if(!is.null(cfh_names)){if(cfh_names!="1") comb<-c(comb,cfh_names)}
   if(length(comb)>1){
     if(length(unique(comb))>1) stop("The same group variable must be used in params, within_grp_var and cfh")
   }
   
  #Substitute actual grouping variable name with 'group'
  if(!is.null(group_name)){
    firstdata<-data[,c(as.character(formula[[2]]),vars)]
    names(firstdata)[which(colnames(firstdata)==group_name)]<-"group"
    firstdata<-firstdata[order(firstdata$group),]
  } else {firstdata<-as.data.frame(data[,c(as.character(formula[[2]]),vars)]); 
     names(firstdata)[1:length(c(formula[[2]],vars))]<-c(as.character(formula[[2]]),vars)}

  #Get model matrix - creates columns for equation use
  if(!is.null(group_name)){ 
    model.temp<-paste("lm(",vars[1],"~as.factor(group),data=firstdata)",sep="")
    cat<-as.data.frame(model.matrix(eval(parse(text=model.temp))))
    if(is.numeric(firstdata$group)) names(cat)<-levels(as.factor(paste("group",firstdata$group,sep="")))
    if(!is.numeric(firstdata$group)|is.factor(firstdata$group)) names(cat)<-paste(group_name,levels(as.factor(firstdata$group)),sep="")
    zz<-as.data.frame(cbind(firstdata[,c(1:length(c(formula[[2]],rhs_equation_vars)))],cat))
    index<-which(is.na(zz),arr.ind=TRUE)
    index<-as.numeric(index[,1])
    if(length(index)>0) zz<-zz[-index,] #remove missing values
  }
  if(is.null(group_name)){ # no group used is equation
    firstdata$group<-1
    zz<-firstdata
  }

   
   
  #Create parameter = values for objective function
  n_prim_parms<-NULL
  for(j in 1:length(names(params))){
    ptnames<-names(params)[j]
    tstart<-eval(parse(text=paste("start$",ptnames,sep="")))
    if(j==1){
      pcnt<-length(tstart)
      n_prim_parms<-c(n_prim_parms,length(tstart))
      ntemp<-length(tstart)
      labL<-paste(ptnames,1:ntemp,"=",sep="")
      labL<-paste(labL[1:ntemp],"parms[",1:ntemp,"]",";",sep="",collapse="")
      parm_labels<-paste(c(labL),collapse="",sep="")
      pcnt<-pcnt+1
    }  
    if(j>1){
      ntemp<-length(tstart)
      n_prim_parms<-c(n_prim_parms,length(tstart))
      tL<-paste(ptnames,1:ntemp,"=",sep="")
      tL<-paste(tL[1:ntemp],"parms[",pcnt:c(pcnt+ntemp-1),"]",";",sep="",collapse="")
      parm_labels<-paste(parm_labels,tL,collapse="",sep="")
      pcnt<-pcnt+ntemp
    }
  }
  
  if(nthetas>0 & is.null(cfh_fixed)){
    labtheta<-paste("theta",1:nthetas,"=",sep="")
    labtheta<-paste(labtheta[1:nthetas],"parms[",pcnt:c(pcnt+nthetas-1),"]",";",sep="") 
    parm_labels<-paste(c(parm_labels,labtheta),collapse="",sep="")
  }
  newformula<-formula
  
 #Manipulate equation include group comparisons
  if(!is.null(group_name)){
    tempdata<-zz
    cat_names<-names(cat)
  # Derive labels for equation variables
  for(j in 1:length(names(params))){
    ptnames<-names(params)[j]
    tstart<-eval(parse(text=paste("start$",ptnames,sep="")))
    pL<-paste(ptnames,1:ngroups,"*","tempdata$",names(tempdata)[length(c(formula[[2]],rhs_equation_vars))+1:c(ncol(tempdata))],sep="")
    if(length(tstart)==1) pL<-pL[1] else pL<-paste(c(paste(pL[1:c(ngroups-1)],"+",sep=""),pL[ngroups]),collapse="")
    newformula[[3]]<-do.call("substitute", list(newformula[[3]], setNames(list(str2lang(pL)), ptnames))) 
  }  
  # Add tempdata to predictors in equation
  for(k in 1:length(rhs_equation_vars)){
    pL<-paste("tempdata$",rhs_equation_vars[k],sep="")
    newformula[[3]]<-do.call("substitute", list(newformula[[3]], setNames(list(str2lang(pL)), rhs_equation_vars[k]))) 
  }
  }

  if(is.null(group_name)){
    tempdata<-zz
    # Derive labels for equation variables
    for(j in 1:length(names(params))){
      ptnames<-names(params)[j]
      tstart<-eval(parse(text=paste("start$",ptnames,sep="")))
      pL<-paste(ptnames,"1","*","tempdata$group",sep="")
      newformula[[3]]<-do.call("substitute", list(newformula[[3]], setNames(list(str2lang(pL)), ptnames))) 
    }
    # Add tempdata to predictors in equation
    for(k in 1:length(rhs_equation_vars)){
      pL<-paste("tempdata$",rhs_equation_vars[k],sep="")
      newformula[[3]]<-do.call("substitute", list(newformula[[3]], setNames(list(str2lang(pL)), rhs_equation_vars[k]))) 
    }
  }
  
###----------------Start to assemble code for fitting -------------------------------###
  tempdata$group<-firstdata$group
  tempdata$theta<-0  
  tempdata$newgroup<-paste(group_name,tempdata$group,sep='')
  if(is.null(cfh)) tempdata$thetagroup<-1
  
  if(is.name(cfh$form[[2]])) tempdata$thetagroup<-as.numeric(as.factor(eval(parse(text=paste("data$",cfh$form[[2]],sep=""))))) #check to make group codes match tempdata
  if(is.numeric(cfh$form[[2]])) tempdata$thetagroup<-1 #check to make group codes match tempdata
  if(is.name(within_grp_var[[2]])) tempdata$sigmagroup<-as.numeric(as.factor(eval(parse(text=paste("data$",within_grp_var[[2]],sep=""))))) #check to make group codes match tempdata
  if(is.numeric(within_grp_var[[2]])) tempdata$sigmagroup<-1 #check to make group codes match tempdata
  
  diff1<-sub("len",as.character(newformula[[2]]),"tempdata$diff1<<-(tempdata$len-tempdata$pred)^2;")#
  sigmastore<-NULL
  LL<-NULL
# No correction for heterogeneity
  if(is.null(cfh)){
     if(nsigmas==1){
       ll_code<-paste(diff1,
                  "ns<-nrow(tempdata);",
                   "sigma2<-sum(tempdata$diff1)/ns;",
                   "sigmastore<<-sigma2;",
                   "LL<--(ns/2)*log(2*pi*sigma2)-(ns/2)",
                   sep="",collapse="")
     }
     if(nsigmas>1){
       ll_code<-paste(diff1,
                   "overall<<-aggregate(tempdata$diff1,list(tempdata$sigmagroup),function(x){c(length(x),sum(x)/length(x))});",
                   "all_comps<-overall$x;",
                   "sigmastore<<-data.frame(overall$Group.1,overall$x[,2]);",
                   "LL<-apply(all_comps,1,function(x){-(x[1]/2)*log(2*pi*x[2])-(x[1]/2)});",
                   "LL<-sum(LL)",sep="",collapse="")
     }
 }

# Correction for heterogeneity
if(!is.null(cfh)){
  thetastore<-NULL
  if(nsigmas==1 & nthetas>0){ #same sigma for all groups,
     if(is.null(cfh_fixed)){ #thetas estimated
       ll_code<-paste(
        diff1,
        "ns<-nrow(tempdata);",
        "if(nthetas==1) tempdata$theta<<-theta1;",
        "if(nthetas>1){",
        "for(i in 1:nthetas) tempdata$theta<<-ifelse(tempdata$thetagroup==i,eval(parse(text=paste('theta',i,sep='',collapse=''))),tempdata$theta);",
        "};",
        "thetastore<<-aggregate(tempdata$theta,list(tempdata$thetagroup),unique);",
        "sigma2<-sum(tempdata$diff1/(ns*abs(tempdata$pred)^(2*tempdata$theta)));",
        "sigmastore<<-sigma2;",
        "LL<--(ns/2)*log(2*pi*sigma2)-sum(tempdata$theta*log(abs(tempdata$pred)))-(ns/2);",sep='',collapse='')
     }
    if(!is.null(cfh_fixed)){ #thetas fixed
      ll_code<-paste(
        diff1,
        "tempdata$theta<<-cfh_fixed[1];",
        "thetastore<<-aggregate(tempdata$theta,list(tempdata$newgroup),unique);", 
        "sigma2<-sum(tempdata$diff1/(nrow(tempdata)*abs(tempdata$pred)^(2*tempdata$theta)));",
        "sigmastore<<-sigma2;",
        "LL<--(nrow(tempdata)/2)*log(2*pi*sigma2)-sum(tempdata$theta*log(abs(tempdata$pred)))-(nrow(tempdata)/2);",sep='',collapse='')
    }
  }
  
  if(nsigmas>1 & nthetas>0){ #different sigma for each groups, 
    if(is.null(cfh_fixed)){ #theta estimated
     sigmastore<-data.frame(group=rep("fill",nsigmas),sigma2=rep("fill",nsigmas))
     ll_code<- paste(
      diff1,
      "if(nthetas==1) tempdata$theta<<-theta1;",
      "if(nthetas>1){",
      "for(i in 1:nthetas) tempdata$theta<<-ifelse(tempdata$thetagroup==i,eval(parse(text=paste('theta',i,sep='',collapse=''))),tempdata$theta);",
      "};",
      "thetastore<<-aggregate(tempdata$theta,list(tempdata$thetagroup),unique);", 
      "tempdata$predtheta<<-abs(tempdata$pred)^(2*tempdata$theta);",
      "lltemp<-NULL;",
      "for(j in 1:nsigmas){",
      "tempo<-tempdata[tempdata$sigmagroup==j,];",
      "n<-nrow(tempo);",
      "sig2<-sum(tempo$diff1/(n*tempo$predtheta));",
      "sigmastore[j,]<<-c(unique(tempo$sigmagroup),sig2);",
      "lltemp<-c(lltemp,-(n/2)*log(2*pi*sig2)-sum(tempo$theta*log(abs(tempo$pred)))-(n/2));};",
      "LL<-sum(lltemp)",sep='',collapse='')
    }
    
    if(!is.null(cfh_fixed)){ #thetas fixed
      sigmastore<-data.frame(group=rep("fill",nsigmas),sigma2=rep("fill",nsigmas))
      ll_code<- paste(
        diff1,
        "if(nthetas==1) tempdata$theta<<-cfh_fixed[1];",
        "if(nthetas>1){",
        "for(i in 1:nthetas) tempdata$theta<<-ifelse(tempdata$thetagroup==i,eval(parse(text=paste('cfh_fixed[',i,']',sep='',collapse=''))),tempdata$theta);",
        "};",
        "thetastore<<-aggregate(tempdata$theta,list(tempdata$thetagroup),unique);", 
        "tempdata$predtheta<<-abs(tempdata$pred)^(2*tempdata$theta);",
        "lltemp<-NULL;",
        "for(j in 1:nsigmas){",
        "tempo<-tempdata[tempdata$sigmagroup==j,];",
        "n<-nrow(tempo);",
        "sig2<-sum(tempo$diff1/(n*tempo$predtheta));",
        "sigmastore[j,]<<-c(unique(tempo$sigmagroup),sig2);",
        "lltemp<-c(lltemp,-(n/2)*log(2*pi*sig2)-sum(tempo$theta*log(abs(tempo$pred)))-(n/2));};",
        "LL<-sum(lltemp)",sep='',collapse='')
    }
  } #nsigmas>1 & nthetas>0
}#correction for heterogeneirty
  
  #Group starting value
    parms<-NULL
    for(j in 1:length(names(params))){
      ptnames<-names(params)[j]
      tstart<-eval(parse(text=paste("start$",ptnames,sep="")))
      parms<-c(parms,tstart)
    }
    if(!is.null(cfh)){
     if(is.null(cfh_fixed)) parms<-c(parms,cfh_value)
    }
    equat<-deparse(newformula[[3]])
  
    
  fit_model<-function(parms){
    eval(parse(text=parm_labels)) # assign parm estimates to equation parms
    tempdata$pred<<-eval(parse(text=equat)) #calculate predicted values
    eval(parse(text=ll_code))
    return(-LL)
  }
  
  output<-try(nlminb(parms,fit_model,control=nlminb.control),silent=TRUE)
  output1<-data.frame()
  class(output1)<-"try-error"
  if(!is(output,"try-error")){
     output1<-try(optim(output$par,fit_model,hessian=TRUE,control=optim.control),silent=TRUE)
  }
  
if(!is(output1,"try-error")){
      varcov<-solve(output1$hessian)
      SE <- sqrt(diag(varcov))
      AIC<--2*(-output1$value)+2*(length(output1$par)+nsigmas) #need to add 1 extra df for each sigma2 
      BIC<--2*(-output1$value)+(length(output1$par)+nsigmas)*log(nrow(tempdata))
      logLik<--output1$value
      conv_info<-data.frame(convergence_code=output1$convergence,message=ifelse(is.null(output1$message),NA,output1$message),
                            function_evaluation=as.numeric(output1$counts[1]),gradient_evaluation=as.numeric(output1$counts[2]))
    compare_df<-length(output1$par)+nsigmas
    df<-nrow(tempdata)-length(output1$par)
    corr_mat<-cov2cor(varcov)
    nvalues<-NULL
    n_model_parms<-length(output1$par)-nthetas #
    nvalues<-n_model_parms+nsigmas+nthetas
    if(!is.null(group_name)) nvalues<-nvalues+ngroups #group par covered by estimatein parms, use ngroups to add n for each group  
    if(is.null(group_name)){if(nsigmas>1){addsigma<-nsigmas;nvalues<-nvalues+addsigma}} #if grp sigmas then will be ns for each grp
    
    nvalues<-nvalues+5#for Loglik,AIC, BIC, Res.df, Res.SE
     group.names<-1
    if(!is.null(group_name)|nsigmas>1|nthetas>1){
      if(!is.null(group_name)) group.names<-paste(grp_name,as.character(levels(as.factor(eval(parse(text=paste("data$",grp_name,sep="")))))),sep="")
      if(nsigmas>1) group.names<-paste(grp_name,as.character(levels(as.factor(eval(parse(text=paste("data$",within_grp_var[[2]],sep="")))))),sep="")
      if(nthetas>1) group.names<-paste(grp_name,as.character(levels(as.factor(eval(parse(text=paste("data$",cfh$form[[2]],sep="")))))),sep="")
    }
  
    name_parms<-NULL
      for(p in 1:length(n_prim_parms)){
        if(p==1){
            if(n_prim_parms[p]==1) name_parms<-names(params)[p] else {
              for(i in 1:n_prim_parms[p]){
                if(i==1)  name_parms<-c(name_parms,paste(names(params)[p],".(Intercept)",sep=""))
                if(i>1)  name_parms<-c(name_parms,paste(names(params)[p],".",group.names[i],sep=""))
             }
            }
        }
        if(p>1){
          if(n_prim_parms[p]==1) name_parms<-c(name_parms,names(params)[p]) else {
            for(i in 1:n_prim_parms[p]){
              if(i==1)  name_parms<-c(name_parms,paste(names(params)[p],".(Intercept)",sep=""))
              if(i>1)  name_parms<-c(name_parms,paste(names(params)[p],".",group.names[i],sep=""))
            }
          }
        }
        
      }
      
    if(nthetas>0){
      if(is.null(cfh_fixed)){
        if(nthetas==1) name_parms<-c(name_parms,"theta") else {
          for(i in 1:nthetas) name_parms<-c(name_parms,paste("theta.",group.names[i],sep=""))}
        colnames(varcov)<-name_parms;rownames(varcov)<-name_parms
        colnames(corr_mat)<-name_parms;rownames(corr_mat)<-name_parms
      }
      if(!is.null(cfh_fixed)){
        colnames(varcov)<-name_parms;rownames(varcov)<-name_parms
        colnames(corr_mat)<-name_parms;rownames(corr_mat)<-name_parms
        if(nthetas==1) name_parms<-c(name_parms,"theta") else {
          for(i in 1:nthetas) name_parms<-c(name_parms,paste("theta.",group.names[i],sep=""))
        }
      }
    }
    if(nthetas<=0){
      colnames(varcov)<-name_parms;rownames(varcov)<-name_parms
      colnames(corr_mat)<-name_parms;rownames(corr_mat)<-name_parms
      
    }
    corr_mat[upper.tri(corr_mat,diag=TRUE)]<-NA
    corr_mat<-round(corr_mat[-1,],5)
    corr_mat<-corr_mat[,-ncol(corr_mat)]
    ns<-NULL
    
    if(!is.null(group_name)|nsigmas>1|nthetas>0){
      dd<-0
      if(nsigmas>1){
        ns<-aggregate(eval(parse(text=sub("len",as.character(newformula[[2]]),"tempdata$len"))),list(as.numeric(as.factor(tempdata$sigmagroup))),length)  #
        names(ns)<-c("group","n")
        dd<-1
      }  
      if(!is.null(group_name)){
        if(dd==0){
          ns<-aggregate(eval(parse(text=sub("len",as.character(newformula[[2]]),"tempdata$len"))),list(as.numeric(as.factor(tempdata$newgroup))),length)  #
          names(ns)<-c("group","n")
          dd<-1
        }
      }
      if(nthetas>0){
        if(dd==0){
          ns<-aggregate(eval(parse(text=sub("len",as.character(newformula[[2]]),"tempdata$len"))),list(as.numeric(as.factor(tempdata$thetagroup))),length)  #
          names(ns)<-c("group","n")
          dd<-1
        }
      }
     
      }

    if(nthetas>0){
        th<-aggregate(tempdata$theta,list(tempdata$thetagroup),unique)
        names(th)<-c("group","theta")
    }
    if(nthetas==0){
        th<-data.frame(group=1,theta=0)
      } 
    n_names<-NULL
    outstats<-NULL
    if(nsigmas>1){
      for(i in 1:nsigmas) name_parms<-c(name_parms,paste("ML.sigma2.",group.names[i],sep=""))
      sigmastore[,2]<-as.numeric(sigmastore[,2])
      names(sigmastore)<-c("group","ML.sigma2")
      sigmas<-sigmastore[,2]
      outstats<-merge(ns,sigmastore,by.x="group",by.y="group",all.x=T,all.y=T)
    }   
    
    if(nsigmas==1){ #there could be group variable for equation so need to create d.f for merging
      name_parms<-c(name_parms,"ML.sigma2")
      sigmas<-as.numeric(sigmastore)
      tsigma<-data.frame(group=1,sigma2=sigmas)
      
      if(!is.null(ns)){
         if(nrow(ns)>1){
           for(kk in 1:c(nrow(ns)-1)){
             tt<-tsigma[1,]
            tsigma<-rbind(tsigma,tt)
             
           } 
           tsigma$group<-c(1:nrow(ns))
         }
      }
      if(is.null(ns)) ns<-data.frame(group=1,n=length(eval(parse(text=sub("len",as.character(newformula[[2]]),"tempdata$len")))))
       outstats<-merge(ns,tsigma,by.x="group",by.y="group",all.x=T,all.y=T)
    }

    if(nrow(outstats)==nrow(th)) outstats<-merge(outstats,th,by.x="group",by.y="group",all.x=T,all.y=T) 
     if(nrow(outstats)!=nrow(th)){
        if(nrow(th)<nrow(outstats)){
          for(kk in 1:c(nrow(outstats)-1)){
              tt<-th[1,]
              th<-rbind(th,tt)
            }
          }
          th$group<-c(1:nrow(outstats))
       outstats<-merge(outstats,th,by.x="group",by.y="group",all.x=T,all.y=T) 
     }
      if(!is.null(grp_name)){
        for(i in 1:nrow(outstats)) n_names<-c(n_names,paste("n.",group.names[i],sep=""))
      }
      if(is.null(grp_name)) n_names<-c(n_names,"n")
    
    name_parms<-c(name_parms,"LogLik","AIC","BIC",paste(n_names,sep=""),"Res. df","Res. SE")
    result <- matrix(NA, ncol = 4, nrow = length(name_parms))
    dimnames(result)[[1]]<-name_parms
    dimnames(result)[[2]]<-c("Value","Std.Error","t-value","Pr(>|t|)")
    
    result[1:nrow(varcov),2]<-round(SE,7)
    ts<-round(output1$par/SE,5)
    pr<-round(2*(1-pt(abs(ts),df)),5) #should not include sigma2
    result[1:nrow(varcov),3]<-ts
    result[1:nrow(varcov),4]<-pr
    
     if(is.null(grp_name)) raw_res<-data.frame(group=1,
                        raw.res=c(eval(parse(text=paste("tempdata$",newformula[[2]],sep="")))-tempdata$pred),
                        fitted=tempdata$pred)
    
    if(!is.null(grp_name)){
      raw_res<-data.frame(group=as.numeric(as.factor(eval(parse(text=paste("data$",grp_name,sep=""))))),
                          raw.res=c(eval(parse(text=paste("tempdata$",newformula[[2]],sep="")))-tempdata$pred),
                          fitted=tempdata$pred)
    }
    tmerge<-merge(raw_res,outstats,by.x="group",by.y="group",all.x=T,all.y=T)
    res_var<-sum((tmerge$raw.res^2)/(df*tmerge$fitted^(2*tmerge$theta)))
    tmerge$stand.res<-tmerge$raw.res/(sqrt(res_var*tmerge$fitted^(2*tmerge$theta)))
    
    if(is.null(cfh)) result[,1]<-round(c(output1$par,sigmas,logLik,AIC,BIC,ns[,2],df,sqrt(res_var)),5)
    if(!is.null(cfh)){
       if(is.null(cfh_fixed)) result[,1]<-round(c(output1$par,sigmas,logLik,AIC,BIC,ns[,2],df,sqrt(res_var)),5)
       if(!is.null(cfh_fixed)) result[,1]<-round(c(output1$par,cfh_fixed,sigmas,logLik,AIC,BIC,ns[,2],df,sqrt(res_var)),5)
    }
    modtype<-"growthlrt_plus"
    fit_tit<-paste("Maximum Likelihood Nonlinear Fit of Model: ",deparse1(formula))
    
    
    nlsout<-list(fitting_method=fit_tit,result,varcov,corr_mat,tmerge[,c("raw.res","stand.res")],tmerge[,c("fitted")],conv_info,compare_df,modtype)
    names(nlsout)<-c("model","results","variance.covariance","correlation","residuals","fitted","convergence","model_comp_df","type")
    return(nlsout)
  }
  if(is(output1,"try-error")) return("Fit Failed.")
  #return(nlsout)
}
