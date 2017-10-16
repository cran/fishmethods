plot.grotagplus <-
function(x,plot.type="meangrowth",Linitial=NULL,
             resid.spec=list(Pearson=T,x='mean.delL'),xlim=NULL,ylim=NULL,
             pch=20,leg.loc=NULL,age.based.growth=NULL,...)
{
    is.in <- function(x, y)!is.na(match(x, y))
    if(is.null(x$datasetnames))
        Ndataset <-1
    else{
        datasetnames <- x$datasetnames
        Ndataset <- length(datasetnames)
    }
    if(is.null(leg.loc))
        leg.loc <- if(plot.type=='traj')'topleft' else 'topright'
    Nrecord <- nrow(x$pred)
    if(plot.type=="meangrowth"){
        pldat <- if(Ndataset==1)list(x$meananngrowth) 
                 else x$meananngrowth
        if(!is.null(age.based.growth)){
            if(!is.list(age.based.growth))
                stop("Argument age.based.growth must be a list")
            if(is.null(names(age.based.growth))){
                if(length(age.based.growth)==Ndataset){
                    if(Ndataset!=1)names(age.based.growth) <- datasetnames
                }
                else stop("If age.based.growth is unnamed it must have",
                          "the same length as x$datasetnames")
            }
            agedat <- lapply(age.based.growth,function(x)
                list(L1=x[-length(x)],delL=diff(x)))
            pldat <- c(pldat,agedat)
            if(is.null(names(age.based.growth))){
                col <- c(1,1)
                lty <- 1:2
                leg <- c('tag','age')
            }
            else if(Ndataset==1){
                col <- 1:length(pldat)
                leg <- c('tag',paste(names(agedat),'age'))
                lty <- c(1,rep(2,length(agedat)))
            }
            
            else{
                extranames <- names(agedat[!is.in(names(agedat),
                                                  datasetnames)])
                col <- c(1:Ndataset,match(names(age.based.growth),
                                          c(datasetnames,extranames)))
                leg <- c(paste(datasetnames,'tag'),paste(names(agedat),'age'))
                lty <- rep(1:2,c(Ndataset,length(agedat)))
            }
        }
        else{
            agedat <- NULL
            col <- 1:Ndataset
            leg <- if(Ndataset>1)datasetnames else NULL
            lty <- rep(1,Ndataset)
        }
        if(is.null(xlim))xlim <- range(unlist(lapply(pldat,function(x)x$L1)))
        if(is.null(ylim))
            ylim <- range(c(0,unlist(lapply(pldat,function(x)x$delL))))
        par(mar=c(2.5,3,1,1)+.1,mgp=c(0,0.5,0),las=1)
        plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='',ylab='')
        for(i in 1:length(pldat))
            lines(pldat[[i]]$L1,pldat[[i]]$delL,col=col[i],lty=lty[i])
        if(length(pldat)>1)legend(leg.loc,leg,lty=lty,col=col,bty='n')
        mtext('Length',1,1.5)
        mtext('Mean annual growth',2,2,las=0)
    }
    else if(plot.type=='traj'){
        if(is.null(Linitial))
            stop("Must provide argument Linitial to plot trajectory")
        nam <- paste('L1=',Linitial,sep='')
        if(is.null(x$traj[[nam]]))
            stop("No trajectory available for Linitial = ",Linitial,
                 "; rerun grotagplus with traj.Linit = ",Linitial)
        xx <- x$traj$T2
        yy <- if(Ndataset==1)list(x$traj[[nam]]) else x$traj[[nam]]
        if(is.null(xlim))xlim <- range(xx)
        if(is.null(ylim))ylim <- range(unlist(yy))
        par(mar=c(2.5,3,1,1)+.1,mgp=c(0,0.5,0),las=1)
        plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='',ylab='')
        for(i in 1:Ndataset)
            for(j in 1:5)lines(xx,yy[[i]][,j],col=i,lty=c(1,2,2,3,3)[j])
        leg <- c('Mean','95% c.i. w/o meas. err.','95% c.i. w. meas. err.')
        lty <- 1:3
        col <- rep(1,3)
        if(Ndataset>1){
            leg <- c(names(yy),leg)
            lty <- c(rep(1,Ndataset),lty)
            col <- c(1:Ndataset,col)
        }
        legend(leg.loc,leg,lty=lty,col=col,bty='n')
        mtext('Time (y)',1,1.5)
        mtext('Length',2,2,las=0)
    }
    else if(plot.type=='resid'){
        yy <- x$pred[[if(resid.spec$Pearson)'Pearson' else 'resid']]
        if(is.in(resid.spec$x,c('L1','delT','mean.delL')))
            xx=x$pred[[resid.spec$x]]
        else stop('Invalid value for resid.spec$x')
        if(is.null(xlim))xlim <- range(xx)
        if(is.null(ylim))ylim <- c(-1,1)*max(abs(unlist(yy)))
        if(is.null(pch))pch <- 20
        par(mar=c(2.5,3,1,1)+.1,mgp=c(0,0.5,0),las=1)
        plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab='',ylab='')
        for(i in 1:Ndataset){
            sel <- (if(Ndataset>1)x$pred$dataset==datasetnames[i]
                    else rep(T,Nrecord))
            points(xx[sel],yy[sel],pch=pch,col=i)
            lines(lowess(xx[sel],yy[sel]),lty=2,col=i)
        }
        xlab <- c(L1='Length at tagging',delT='Time at liberty (y)',
                  mean.delL='Expected length increment')[resid.spec$x]
        ylab <- ifelse(resid.spec$Pearson,'Pearson residual','Residual')
        mtext(xlab,1,1.5)
        mtext(ylab,2,2,las=0)
        abline(h=0)
        if(resid.spec$Pearson)abline(h=c(-2,2),lty=3)
        if(Ndataset>1)legend(leg.loc,datasetnames,pch=pch,col=1:Ndataset,
                             bty='n')
    }
    else stop("Invalid value for argument plot.type")
    invisible()
}

