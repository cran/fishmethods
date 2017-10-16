print.grotagplus <-
function(x,precision=c(est='sig3',stats='dec1',cor='dec2'),...)
{
    is.in <- function(x, y)!is.na(match(x, y))
    reduce <- function(comp,prec){
        if(substring(prec,1,3)=='sig')
            signif(comp,as.numeric(substring(prec,4)))
        else  if(substring(prec,1,3)=='dec')
            round(comp,as.numeric(substring(prec,4)))
        else stop('Invalid value for precision: "',prec,'"')
    }
    conv <- c(parest='est',parfix='est',stats='stats',
              correlations='cor',Linf.k='est')
    pr.comp <- c('model','parest',"parfix",'stats','correlations',"Linf.k")
    pr.comp <- pr.comp[is.in(pr.comp,names(x))]
    for(nam in pr.comp[-1])
        x[[nam]] <- reduce(x[[nam]],precision[conv[nam]])
    print.default(x[pr.comp])
}
