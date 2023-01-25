grotagplus <-
function(tagdata,dataID=NULL,alpha,beta=NULL,
             model=list(mean='Francis',var='linear',seas='sinusoid'),
             design,stvalue,upper,lower,fixvalue=NULL,traj.Linit=c(alpha,beta),
             control=list(maxit=10000),debug=FALSE)
{
    is.in <- function(x,y)!is.na(match(x,y))
    getparlist <- function(p)
        ## Convert the parameter vector to a list of vectors so that, e.g.,
        ## pars$galpha[i] is the value of galpha appropriate for the ith data 
        ## record, adding fixed and default parameter values as necessary
    {
        pars <- list() 
        i1 <- 1
        for(nam in est.parnames){
            i2 <- i1+Nvalue[[nam]]-1
            pars[[nam]] <- p[i1:i2][parindex[[nam]]]
            i1 <- i1+Nvalue[[nam]]
        }
        for(nam in fix.parnames){
            pars[[nam]] <- (if(!is.list(fixvalue[[nam]]))
                                fixvalue[[nam]][parindex[[nam]]]
                            else
                                fixvalue[[nam]]$value[parindex[[nam]]])
        }
        for(nam in names(default.parvals)) ## add default values if necessary
            if(!is.in(nam,names(pars)))pars[[nam]] <- default.parvals[[nam]]
        pars
    }
    phi <- function(T,pars){
        if(model$seas=='sinusoid')pars$u*(sin(2*pi*(T-pars$w)))/(2*pi)
        else if(model$seas=='switched'){
            z <- (T - pars$w - 0.5/(pars$u + 1))%%1
            ifelse(z<=(pars$u/(pars$u+1)),0.5*(pars$u/(pars$u+1))-z,
                   pars$u*(z-(pars$u+0.5)/(pars$u+1)))
        }
        else 0
    }
    getmn <- function(pars,tdiff,L1){ # get expected length increment
        if(model$mean=='Francis' | model$mean=='asymptotic'){
            Linf <- (beta * pars$galpha - alpha * pars$gbeta)/(pars$galpha - pars$gbeta)
            vonbk <- -log(1 + (pars$galpha - pars$gbeta) /(alpha - beta))
            if(model$mean=='Francis')
                mu <- (Linf - L1) * (1 - exp(-vonbk*tdiff))
            else{ # model = 'asymptotic' 
                mu <- rep(0,length(L1))
                OK <- pars$Lstar<Linf
                if(any(OK)){
                    Lstar1 <- Linf-(Linf-pars$Lstar)*exp(vonbk*tdiff)
                    L1group <- ifelse(L1<=Lstar1,'lo',
                               ifelse(L1>=pars$Lstar,'hi','med'))
                    sel <- OK & L1group=='lo'
                    if(any(sel))
                        mu[sel] <- (Linf[sel]- L1[sel]) *
                            (1 - exp(-vonbk[sel]*tdiff[sel]))
                    sel <- OK & L1group=='med'
                    if(any(sel))
                        mu[sel] <- pars$Lstar[sel]-L1[sel] +
                            (Linf[sel]-pars$Lstar[sel])*
                            log(1+vonbk[sel]*tdiff[sel]-
                                log((Linf[sel]-L1[sel])/
                                    (Linf[sel]-pars$Lstar[sel])))
                    sel <- OK & L1group=='hi'
                    if(any(sel))
                        mu[sel] <- (Linf[sel]-pars$Lstar[sel])*
                            log(1+vonbk[sel]*tdiff[sel]*
                                exp((pars$Lstar[sel]-L1[sel])/
                                    (Linf[sel]-pars$Lstar[sel])))
                }
                else{
                    sel <- !OK & L1<Linf
                    if(any(sel))
                        mu[sel] <- (Linf[sel] - L1[sel]) *
                            (1 - exp(-vonbk[sel]*tdiff[sel]))
                }
            }
        }
        else if(model$mean=='Schnute' | model$mean=='Schnute.aeq0'){
            lama <- alpha+pars$galpha
            beq0 <- pars$b==0
            a <- c <- mu <- rep(0,length(beq0))
            if(model$mean=='Schnute'){
                lamb <- beta+pars$gbeta
                a[!beq0] <- log((beta^pars$b[!beq0]-alpha^pars$b[!beq0])/
                                (lamb[!beq0]^pars$b[!beq0]-
                                 lama[!beq0]^pars$b[!beq0]))
                a[beq0] <- log((log(beta/alpha))/(log(lamb[beq0]/lama[beq0])))
                c[!beq0] <- ((beta*lama[!beq0])^pars$b[!beq0] -
                             (alpha*lamb[!beq0])^pars$b[!beq0])/
                    (lama[!beq0]^pars$b[!beq0]-alpha^pars$b[!beq0]+
                     beta^pars$b[!beq0]-lamb[!beq0]^pars$b[!beq0])
                c[beq0] <- (log(beta)*log(lama[beq0])-log(alpha)*log(lamb[beq0]))/
                    (log(lama[beq0]*beta)-log(lama[beq0]*alpha))
            }
            aeq0 <- a==0
            sel <- !aeq0 & !beq0
            if(any(sel))
                mu[sel] <- -L1[sel]+((L1[sel]^pars$b[sel])*
                                     exp(-a[sel]*tdiff[sel]) +
                                     c[sel]*(1-exp(-a[sel]*tdiff[sel])))^
                    (1/pars$b[sel])
            sel <- aeq0 & !beq0
            if(any(sel))
                mu[sel] <- -L1[sel]+((L1[sel]^pars$b[sel])+
                                     (lama[sel]^pars$b[sel]-
                                      alpha^pars$b[sel])*tdiff[sel])^
                    (1/pars$b[sel])
            sel <- !aeq0 & beq0
            if(any(sel))
                mu[sel] <- -L1[sel]+(L1[sel]^
                                     exp(-a[sel]*tdiff[sel]))*
                    exp(c[sel]*(1-exp(-a[sel]*tdiff[sel])))
            sel <- aeq0 & beq0
            if(any(sel))
                mu[sel] <- -L1[sel]+L1[sel]*
                    (lama[sel]/alpha)^tdiff[sel]
        }
        as.vector(mu)
    }
    getsd <- function(mn,pars){ # get s.d. of the length increment
        mn <- ifelse(mn<0,1,mn)
        if(model$var=='linear')
            sigma <- pars$nu * mn
        else if(model$var=='capped'){
            sigma <- pars$nu * mn
            sigma[sigma>pars$t] <- pars$t
        }
        else if(model$var=='exponential')
            sigma <- pars$t*(1-exp(-pars$nu * mn))
        else if(model$var=='power')
            sigma <- pars$nu * (mn^pars$t)
        else if(model$var=='least-squares')
            sigma <- 0
        sdev <- sqrt(sigma^2+ pars$s^2)
        as.vector(sdev)
    }
    ## Check alpha and beta
    if(!is.null(beta)){
        if(alpha >= beta)
            stop("Error: parameter alpha must be smaller than beta")
    }
    else if(model$mean!='Schnute.aeq0')
        stop ('argument "beta" is missing')
    ## Check model
    validmodels <- list(mean=c('Francis','Schnute','Schnute.aeq0','asymptotic'),
                        var=c('linear','capped','exponential','power',
                              'least-squares'),
                        seas=c('sinusoid','switched','none'))
    for(nam in names(validmodels)){
        if(!is.in(model[[nam]],validmodels[[nam]]))
            stop('Invalid value for model$',nam,': ',model[[nam]])
    }
    ## Define parnames for selected model
    model.parnames <- c(list(
        Francis=c('galpha','gbeta'),
        Schnute=c('galpha','gbeta','b'),
        Schnute.aeq0=c('galpha','b'),
        asymptotic=c('galpha','gbeta','Lstar'))[[model$mean]],
        if(model$seas!='none')c('u','w'),
        if(model$var=='least-squares')'s'
        else if(model$var=='linear')c('nu','s','m') else c('nu','t','s','m'),
        'p')
    if(debug)cat('model.parnames',model.parnames,'\n')
    ## Parameters which have default values
    default.parvals <- list(nu=0,m=0,p=0)
    ## Check tagdata
    compnames <- c('L1','L2','T1','T2')
    missing <- compnames[!is.in(compnames,names(tagdata))]
    if(length(missing)>0)
        stop("Following component(s)","missing from tagdata:",
             paste(missing,collapse=', '))
    if(!is.null(dataID)){
        if(is.in(dataID,names(tagdata))){
            dataset <- tagdata[[dataID]]
            datasetnames <- as.character(sort(unique(dataset)))
            Ndataset <- length(datasetnames)
        }
        else stop(dataID," is not a component of tagdata")
    }
    else{
        datasetnames <- NULL
        dataset <- rep(1,nrow(tagdata))
        Ndataset <- 1
    }
    ## Drop data records with missing fields and add delL, delT
    tagdata <- tagdata[!apply(is.na(tagdata[,c(compnames,dataID)]),1,any),]
    tagdata$delL <- tagdata$L2-tagdata$L1
    tagdata$delT <- tagdata$T2-tagdata$T1
    Nrecord <- nrow(tagdata)
    ## Check design is valid
    missing <- model.parnames[!is.in(model.parnames,names(design))]
    if(length(missing)>0)
        stop("Following component(s) missing from design:",
             paste(missing,collapse=', '))
    invalid <- names(design)[!is.in(names(design),model.parnames)]
    if(length(invalid)>0)
        stop("Following component(s) of design is/are not for valid",
             "parameter(s):",paste(invalid,collapse=', '))
    Nvalue <- list() ## number of values to be estimated for each parameter
    for(nam in names(design)){
        if(length(design[[nam]])==1){
            if(!is.in(design[[nam]],c(0,1)))
                stop("Value of design$",nam," is invalid",sep='')
            if(design[[nam]]==1){
                Nvalue[[nam]] <- 1
                design[[nam]] <- rep(1,Ndataset)
            }
            else Nvalue[[nam]] <- 0
        }
        else{
            missing <- datasetnames[!is.in(datasetnames,unlist(design[[nam]]))]
            if(length(missing)>0)
                stop("Following dataset(s) missing from design$",nam,":",
                     paste(missing,collapse=', '),sep='')
            uvals <- unique(unlist(design[[nam]]))
            invalid <- uvals[!is.in(uvals,datasetnames)]
            if(length(invalid>0))
                stop("Invalid dataset value(s) in design$",nam,":",
                     paste(invalid,collapse=', '),sep='')
            nval <- table(unlist(design[[nam]]))
            if(any(nval>1))
                stop("Following dataset(s) repeated in design$",nam,":",
                     paste(names(nval)[nval>1],collapse=', '),sep='')
            Nvalue[[nam]] <- length(design[[nam]])
        }
    }
    est.parnames <- names(Nvalue)[unlist(lapply(Nvalue,function(x)x>0))]
    if(debug)cat('Nvalue',unlist(Nvalue),'\n')
    ## Check structure of stvalue, upper, lower
    calc.stvalue <- list(galpha=F,gbeta=F) # Which starting values to calculate
    for(arg.name in c('stvalue','upper','lower')){
        this.arg <- eval(parse(text=arg.name))
        extra <- names(this.arg)[!is.in(names(this.arg),est.parnames)]
        if(length(extra)>0)
            warning("Following unestimated component(s) of ",arg.name,
                    " ignored: ",paste(extra,collapse=', '))
        for(nam in est.parnames){
            this.comp <- this.arg[[nam]]
            if(length(this.comp)==0){
                if(nam=='galpha' | nam=='gbeta'){
                    if(arg.name=='stvalue')
                        calc.stvalue[[nam]] <- T
                    else{
                        if(!calc.stvalue[[nam]])
                            stop("Component ",nam," missing from ",arg.name)
                    }
                }
                else stop("Component ",nam," missing from ",arg.name)
            }
            else if(length(this.comp)==1)
                this.arg[[nam]] <- rep(this.comp,Nvalue[[nam]])
            else if(length(this.comp)!=Nvalue[[nam]])
                stop("Invalid length for ",arg.name,'$',nam,sep='')
        }
    }
    ## Check that upper and lower bracket stvalue
    for(nam in names(stvalue)){
        if(!all(lower[[nam]]<upper[[nam]]))
            stop('lower not less than upper for parameter ',nam)
        if(!all(lower[[nam]]<=stvalue[[nam]]))
            stop('lower not less than or equal to stvalue for parameter ',nam)
        if(!all(upper[[nam]]>=stvalue[[nam]]))
            stop('upper not greater than or equal to stvalue for parameter ',nam)
    }
    ## Check format & content of fixvalue and add fixed parameters to list Nvalue
    fix.parnames <- names(fixvalue)
    invalid <- fix.parnames[!is.in(fix.parnames,model.parnames) |
                            is.in(fix.parnames,est.parnames)]
    if(length(invalid)>0)
        stop("Following component(s) should not be in fixvalue:",
             paste(invalid,collapse=', '))
    missing <- model.parnames[!is.in(model.parnames,c(est.parnames,fix.parnames,
                                                      names(default.parvals)))]
    if(length(missing)>0)
        stop("Following parameter(s) are not estimated and have no\n",
             " default value so they must be included in fixvalue:",
             paste(missing,collapse=', '))
    for(nam in fix.parnames){
        if(!is.list(fixvalue[[nam]])){
            if(!is.numeric(fixvalue[[nam]]) | length(fixvalue[[nam]])>1)
                stop(paste("fixvalue",nam,sep='$'),"is invalid:",
                     "must be a single number")
            Nvalue[[nam]] <- 1
        }
        else{
            missing <- c('design','value')[!is.in(c('design','value'),
                                                  names(fixvalue[[nam]]))]
            if(length(missing)>0)
                stop("Following component(s) missing from list",
                     paste("fixvalue",nam,sep='$'))
            missing <- datasetnames[!is.in(datasetnames,
                                           unlist(fixvalue[[nam]]$design))]
            if(length(missing)>0)
                stop("Following dataset(s) missing from fixvalue$",nam,
                     "$design: ",paste(missing,collapse=', '),sep='')
            uvals <- unique(unlist(fixvalue[[nam]]$design))
            invalid <- uvals[!is.in(uvals,datasetnames)]
            if(length(invalid>0))
                stop("Invalid dataset value(s) in fixvalue$",nam,":",
                     "$design: ",paste(invalid,collapse=', '),sep='')
            nval <- table(unlist(design[[nam]]))
            if(any(nval>1))
                stop("Following dataset(s) repeated in fixvalue$",nam,":",
                     "$design: ",paste(names(nval)[nval>1],collapse=', '),sep='')
            Nvalue[[nam]] <- length(fixvalue[[nam]]$value)
        }
    }
    ## Adjust calc.stvalue if necessary
    for(nam in c('galpha','gbeta'))
        if(is.in(nam,fix.parnames))calc.stvalue[[nam]] <- F
    ## Make parindex, linking data records to parameter values
    ## (e.g., parindex$m[i]==j means that the likelihood calculation for
    ## the ith data record uses the jth value of parameter m)
    ##
    if(debug)cat('est.parnames:',est.parnames,'\n')
    parindex <- list()
    for(nam in est.parnames){## First incude all estimated parameters
        if(Nvalue[[nam]]==1)parindex[[nam]] <- rep(1,Nrecord)
        else{
            parindex[[nam]] <- rep(1,Nrecord)
            for(i in 2:(Nvalue[[nam]]))
                parindex[[nam]][is.in(dataset,design[[nam]][[i]])] <- i
        }
    }
    for(nam in fix.parnames){ ## Then add any fixed values
        if(Nvalue[[nam]]==1)parindex[[nam]] <- rep(1,Nrecord)
        else{
            parindex[[nam]] <- rep(1,Nrecord)
            for(i in 2:(Nvalue[[nam]]))
                parindex[[nam]][is.in(dataset,fixvalue[[nam]]$design[[i]])] <- i
        }
    }
    ## Calculate galpha, gbeta starting values and bounds, if needed
    if(any(unlist(calc.stvalue))){
        for(nam in names(calc.stvalue)[unlist(calc.stvalue)]){
            Lref <- ifelse(nam=='galpha',alpha,beta)
            sel <- abs(tagdata$L1/Lref-1)<0.2 & (tagdata$delT >
                                                 median(tagdata$delT))
            stvalue[[nam]] <- as.vector(tapply(tagdata$delL[sel]/
                                               tagdata$delT[sel],
                                               dataset[sel],mean))
            if(!is.in(nam,names(lower)))
                lower[[nam]] <- 0.5*stvalue[[nam]]
            if(!is.in(nam,names(upper)))
                upper[[nam]] <- 1.5*stvalue[[nam]]
            if(debug)cat(nam,sum(sel),stvalue[[nam]],lower[[nam]],
                         upper[[nam]],'\n')
        }
    }
    startpar <- unlist(stvalue[est.parnames])
    all.pars.fixed <- length(startpar)==0
    lower2 <- unlist(lower[est.parnames])
    upper2 <- unlist(upper[est.parnames])
    if(debug)print(rbind(start=startpar,lo=lower2,hi=upper2))
    ## Negative log-likelihood function
    NLL<-function(p){
        pars <- getparlist(p)
        ##        print(unlapply(pars,mean))
        R <-diff(range(tagdata$delL))
        tdiff <- tagdata$delT+phi(tagdata$T2,pars)-phi(tagdata$T1,pars)
        mu <- getmn(pars,tdiff,tagdata$L1)
        sdev <- getsd(mu,pars)
        NLL <- -1*sum(log((1-pars$p) * dnorm(tagdata$delL,mu+pars$m,sdev) +
                          pars$p/R))
        if(is.infinite(NLL)){
            cat('Log-likelihood is not finite with following pars\n')
            print(p)
        }
        NLL
    }
    if(debug){
        pars <- getparlist(startpar)
        tdiff <- tagdata$delT+phi(tagdata$T2,pars)-phi(tagdata$T1,pars)
        mu <- getmn(pars,tdiff,tagdata$L1)
        sdev <- getsd(mu,pars)
        cat('mu and sdev for startpar\n')
        print(rbind(mu=summary(mu),sdev=summary(sdev)))
        cat("NLL of startpar",NLL(startpar),'\n')
    }
    if(!all.pars.fixed){
        mod1<-try(optim(par=startpar,fn=NLL,lower=lower2,upper=upper2,
                        method="L-BFGS-B",hessian=TRUE, control=control),
                  silent=TRUE)
    }
    else
        mod1 <- list(value=NLL(startpar),convergence=0)
       classtype<-class(mod1)
    if(classtype=="try-error"){
        cat('Fit failed\n')
        stop(mod1$message)
    }
    else{
        if(mod1$convergence!=0)cat('Warning:',mod1$message,'\n')
        out <- list()
        parnames.out <- model.parnames[is.in(model.parnames,est.parnames)]
        ## Make matrix of parameter estimates and their s.e.s
        if(!all.pars.fixed){
            covmat <-  solve(mod1$hessian)
            SE<-sqrt(diag(covmat))
            cormat <- cov2cor(covmat)
            dimnames(cormat) <- list(names(startpar),names(startpar))
            rownam <- (if(Ndataset==1)c('est','se')
                       else c(paste(datasetnames,'.est',sep=''),
                              paste(datasetnames,'.se',sep='')))
            parest <- matrix(0,2*Ndataset,length(est.parnames),
                             dimnames=list(rownam,parnames.out))
            i1 <- 1
            for(nam in est.parnames){
                if(Nvalue[[nam]]==1)parest[,nam] <-
                                       rep(c(mod1$par[i1],SE[i1]),rep(Ndataset,2))
                else{
                    for(j in 1:Nvalue[[nam]]){
                        parest[paste(design[[nam]][[j]],'.est',sep=''),nam] <-
                            mod1$par[i1+j-1]
                        parest[paste(design[[nam]][[j]],'.se',sep=''),nam] <-
                            SE[i1+j-1]
                    }
                }
                i1 <- i1+Nvalue[[nam]]
            }
            out$parest <- parest
        }
        ## Make matrix of fixed parameters (if there are any)
        rownam <- if(Ndataset==1)'' else datasetnames
        if(length(fix.parnames)>0){
            parfix <- matrix(0,Ndataset,length(fix.parnames),
                             dimnames=list(rownam,fix.parnames))
            for(nam in fix.parnames){
                if(Nvalue[[nam]]==1)parfix[,nam] <- fixvalue[[nam]]
                else{
                    for(j in 1:Nvalue[[nam]]){
                        parfix[fixvalue[[nam]]$design[[j]],nam] <-
                            fixvalue[[nam]]$value[j]
                    }
                }
            }
            out$parfix <- parfix
        }
        ## Fit statistics etc
        if(!all.pars.fixed){            
            out$correlations <- cormat
            out$stats <- c(negloglikl=mod1$value,
                           AIC=2*mod1$value+2*length(mod1$par))
        }
        else
            out$stats <- round(c(negloglikl=mod1$value),1)
        out$model <- unlist(model)
        if(!is.null(datasetnames))out$datasetnames <- datasetnames
        ## Make pred - dataframe including stuff needed for residual plots etc
        pars <- getparlist(mod1$par)
        tdiff <- tagdata$delT+phi(tagdata$T2,pars)-phi(tagdata$T1,pars)
        mean.delL <- getmn(pars,tdiff,tagdata$L1)
        sd.delL <- getsd(mean.delL,pars)
        out$pred <- cbind(tagdata[,c('L1','delT')],mean.delL=mean.delL,
                          sd.delL=sd.delL,resid=tagdata$delL-mean.delL)
        out$pred$Pearson <- out$pred$resid/sd.delL
        if(Ndataset>1)out$pred$dataset <- dataset
        ## Calculate Linf,vonbk (where possible)
        if(model$mean!='Schnute.aeq0' & model$mean!='Schnute'){
            galpha <- (if(is.in('galpha',est.parnames))parest[1:Ndataset,'galpha']
                       else parfix[,'galpha'])
            gbeta <- (if(is.in('gbeta',est.parnames))parest[1:Ndataset,'gbeta']
                      else parfix[,'gbeta'])
            out$Linf.k <- cbind(Linf=(beta * galpha - alpha * gbeta)/
                                    (galpha - gbeta),
                                vonbk=-log(1 + (galpha - gbeta)/(alpha - beta)))
            rownames(out$Linf.k) <- if(Ndataset>1)datasetnames else ''
        }
        ## Data for plot of mean annual growth vs L1
        out$meananngrowth <- list()
        nam <- c('L1','delL')
        for(i in 1:Ndataset){
            sel <- if(Ndataset==1)rep(T,Nrecord) else dataset==datasetnames[i]
            L1 <- seq(min(tagdata$L1[sel]),max(tagdata$L1[sel]),length=200)
            if(Ndataset==1){
                pars2 <- lapply(pars,function(x)rep(x[1],length(L1)))
                out$meananngrowth <-
                    as.data.frame(list(L1=L1,
                                       delL=getmn(pars2,rep(1,length(L1)),L1)))
            }
            else{
                pars2 <- lapply(pars,
                                function(x)x[rep(match(datasetnames[i],dataset),
                                                 length(L1))])
                out$meananngrowth[[datasetnames[i]]] <-
                    as.data.frame(list(L1=L1,
                                       delL=getmn(pars2,rep(1,length(L1)),L1)))
            }
        }
        ## 1-year growth trajectories
        T2 <- seq(0,1,1/24)
        traj <- list(T2=T2)
        for(Linit in traj.Linit){
            nam <- paste('L1=',Linit,sep='')
            traj[[nam]] <- list()
            L1 <- rep(Linit,length(T2))
            for(i in 1:Ndataset){
                if(Ndataset==1)
                    pars2 <- lapply(pars,function(x)rep(x[1],length(T2)))
                else
                    pars2 <- lapply(pars,
                                    function(x)x[rep(match(datasetnames[i],
                                                           dataset),length(T2))])
                tdiff <- T2+phi(T2,pars2)-phi(0,pars2)
                mn <- getmn(pars2,tdiff,L1)
                sd.obs <- getsd(mn,pars2)
                sd.gro <- sqrt(sd.obs^2-pars2$s^2)
                mat <- as.data.frame(list(L2=L1+mn,gro.lo=L1+mn-2*sd.gro,
                                          gro.hi=L1+mn+2*sd.gro,
                                          obs.lo=L1+mn-2*sd.obs,
                                          obs.hi=L1+mn+2*sd.obs))
                if(Ndataset==1)traj[[nam]] <- mat
                else traj[[nam]][[datasetnames[i]]] <- mat
            }
        }
        out$traj <- traj
        if(debug)out$mod1 <- mod1
        class(out) <- "grotagplus"
        return(out)
    }
}
