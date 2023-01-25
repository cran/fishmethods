\name{grotagplus}
\alias{grotagplus}
\title{
Flexible maximum likelihood estimation of growth from multiple tagging datasets.
}
\description{
This is an extension of fishmethods function grotag to allow a wider
variety of growth models and also the simultaneous analysis
of multiple tagging datasets with parameter sharing between
datasets (see Details).

As in grotag, the data are fitted using a constrained maximum
likelihood optimization performed by optim using the "L-BFGS-B" method.
Estimated parameters can include galpha, gbeta (mean annual growth at reference
lengths alpha and beta); b (a curvature parameter for the Schnute
models); Lstar (a transitional length for the asymptotic model); m, s (mean and
s.d. of the measurement error for length increment); nu, t (growth
variability); p (outlier probability); u, w (magnitude and phase of
seasonal growth).
}
\usage{
grotagplus(tagdata, dataID=NULL,alpha, beta = NULL,
 model=list(mean="Francis",var="linear",seas="sinusoid"),
 design, stvalue, upper, lower,fixvalue=NULL,
 traj.Linit=c(alpha,beta),control = list(maxit = 10000), debug = FALSE)
}
\arguments{
   \item{tagdata}{Dataframe with components L1, L2 (lengths at release
       and recovery of tagged fish), T1, T2 (julian times (y) at
       release and recovery), and (optionally), a numeric or character
       vector (named by argument \code{dataID}) identifying which
       dataset each data record belongs to (with n datasets this must
       include n unique values). Other components are ignored, as are any records with missing
       values in the required components.}
  \item{dataID}{Name of optional component of tagdata identifying separate
       datasets within tagdata.  The default \code{dataID}=NULL means
       there is no such component (so there is only one dataset).}
  \item{alpha}{Numeric value giving an arbitrary length alpha.}
  \item{beta}{Numeric value giving an arbitrary length beta
       (must have \code{beta} > \code{alpha}).}
  \item{model}{List with components mean, var, seas, specifying which
       model equations to use for the mean (or expected) growth,
       individual variability in growth, and seasonal variation in
       growth (see Details for valid values). The default is that of
       model 4 in Francis (1988).}
  \item{design}{List specifying the design of the estimation: which
       parameters are estimated, and whether multiple values are estimated.
       There should be one component for each parameter of the model
       specified by \code{model}.  Each component must be
       either 0 (not estimated), 1 (same parameter value estimated for all data),
       or, when there are multiple datasets, a list in which each
       component is a sub-vector of unique(tagdata[[dataID]]) and all
       members of unique(tagdata[[dataID]]) occur in one and only one
       component of the list (e.g., galpha=list("Area2",c("Area1", "Area3") )
       means that two values of galpha are to be estimated: one applying to
       the dataset Area2, and the other to datasets Area1 and Area3).}
  \item{stvalue}{List containing starting values of estimated parameters,
       used as input in the nonlinear estimation
       (function \emph{optim}) routine. There should be one component
       for each estimated parameter (except, optionally, galpha and gbeta).
       Each component should be either a single number or a vector whose
       length is the number of separate values of that parameter
       (as specified in \code{design}). In the latter case,
       the order of the parameter values should correspond to that in
       \code{design} (e.g., if design$galpha is as above and
       \code{stvalue}$galpha=c(10,15) then 10 will apply to Area2 and
       15 to Area1 & Area3).  If galpha or gbeta are omitted from stvalue then
       their starting values are calculated from the data.}
   \item{lower}{Lists containing lower limits for each parameter,
       with structure as for \code{stvalue}. galpha and/or gbeta may be omitted
       if they don"t appear in \code{stvalue}.}
   \item{upper}{Lists containing upper limits for each parameter,
       with structure as for \code{stvalue}. galpha and/or gbeta may be omitted
       if they don"t appear in \code{stvalue}.}
  \item{fixvalue}{Optional list containing fixed values for parameters that
       are needed (according to \code{model}) but are not being
       estimated (according to \code{design}) and do not have default
       values (the only default parameter values are nu = 0, m = 0, p = 0). 
       The list should have one named component for each fixed
       parameter.  Usually, each component will be a single number. See
       example below for the required format when a fixed parameter
       takes  different values for different datasets.}
  \item{traj.Linit}{Vector of initial length(s) for output growth trajectories.
       Default is c(alpha,beta).}
  \item{control}{Additional controls passed to the optimization function
       \emph{optim}.}
  \item{debug}{output debugging information.} 
}
\details{
 Valid values of model$mean are
  \code{"Francis"} as in Francis (1988).
  \code{"Schnute"} as in Francis (1995).
  \code{"Schnute.aeq0"} special case of Schnute - see equns (5.3), (5.4)
    of Francis (1995).
  \code{"asymptotic"} as in Cranfield et al. (1996).
    
  Valid values of model$var are
  \code{"linear"} as used in the example in Francis(1988) - see equn
    (5).
  \code{"capped"} as in equn (6) of Francis(1988).
  \code{"exponential"} as in equn (7) of Francis(1988).
  
  \code{"asymptotic"} as in equn (8) of Francis(1988).
  \code{"least-squares"} ignore individual variability and fit data by
    least-squares, as in Model 1 of Francis(1988).  

  Valid values of model$seas are
  \code{"sinusoid"} as in model 4 of Francis(1988).
  \code{"switched"} as in Francis & Winstanley (1989).
  \code{"none"} as in all but model 4 of Francis(1988).
  
  The option of multiple data sets with parameter sharing is intended to
  allow for the situation where we wish to estimate different mean
  growth for two or more datasets but can reasonably assume that
  other parameters (e.g., for growth variability, measurement error,
  outlier contamination) are the same for all datasets.
  This should produces stronger estimates of these other parameters.
  For example, Francis & Francis (1992) allow growth to differ by sex, and
  in Francis & Winstanley (1989) it differs by stock and/or habitat.

  \code{grotagplus} may fail if parameter starting values are too distant from
  their true value, or if parameter bounds are too wide.  Try changing
  these values.  Sometimes reasonable starting values can be found by
  fitting the model with other parameters fixed at plausible values. 
}
\value{
  \item{parest}{Parameter estimates and their s.e.s.}
  \item{parfix}{Parameter values, if any, fixed by user.}
  \item{correlations}{Correlations between parameter estimates. When
    there are multiple estimates of a parameter these are numbered
    by their ordering in argument \code{design}, so in example given
    above galpha1 would apply to Area1, and galpha2 to Area2 and Area3.}
  \item{stats}{Negative log-likelihood and AIC statistic.}
  \item{model}{The three components of the grotagplus argument model.}
  \item{datasetnames}{The dataset names, if there are multiple datasets.}
  \item{pred}{Dataframe of various predicted quantities need for
    residual plots - one row per data record.}
  \item{Linf.k}{Values of parameters Linf and k as calculated between
    equations (1) and (2) of Francis (1988) (but not possible for the
    Schnute model).  These are provided for computational convenience
    only; they are not comparable with Linf and k estimated from
    age-length data. Comparisons of growth estimates from tagging
    and age-length data are better done using output \code{meananngrowth}.}
  \item{meananngrowth}{Data for plot of mean annual growth vs length,
    as in Fig. 8 of Francis and Francis (1992).}
  \item{traj}{Data for plots of growth trajectories like Fig. 2 of
    Francis (1988).}
}
\references{ \code{1} Francis, R.I.C.C., 1988. Maximum likelihood estimation of
  growth and growth variability from tagging data.
  New Zealand Journal of Marine and Freshwater Research, 22, p.42-51.

 \code{2} Cranfield, H.J., Michael, K.P., and Francis, R.I.C.C. 1996.
  Growth rates of five species of subtidal clam on a beach in the South
  Island, New Zealand.  Marine and Freshwater Research 47: 773-784.
  
 \code{3} Francis, R.I.C.C. 1995.  An alternative mark-recapture analogue
 of Schnute"s growth model.  Fisheries Research 23: 95-111.

 \code{4} Francis, R.I.C.C. and Winstanley, R.H. 1989. Differences in
  growth rates between habitats of southeast Australian snapper
  (Chrysophrys auratus). Australian Journal of Marine & Freshwater
  Research 40: 703-710.
  
 \code{5} Francis, M.P. and Francis, R.I.C.C. 1992.  Growth rate
  estimates for New Zealand rig (Mustelus lenticulatus).  Australian
  Journal of Marine and Freshwater Research 43: 1157-1176. 
}
\author{
Chris Francis \email{chrisfrancis341@gmail.com}

Marco Kienzle \email{Marco.Kienzle@gmail.com}

Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\seealso{\code{\link{plot.grotagplus}} \code{\link{print.grotagplus}}}
\examples{
#Model 4 of Francis (1988)
data(bonito)
grotagplus(bonito,alpha=35,beta=55,
               design=list(galpha=1,gbeta=1,s=1,nu=1,m=1,p=1,u=1,w=1),
               stvalue=list(s=0.81,nu=0.3,m=0,p=0.01,u=0.5,w=0.5),
               upper=list(s=3,nu=1,m=2,p=0.1,u=1,w=1),
               lower=list(s=0.1,nu=0.1,m=-2,p=0,u=0,w=0))

#Model 1 of Francis (1988), using least-squares fit
grotagplus(bonito,alpha=35,beta=55,
               model=list(mean="Francis",var="least-squares",seas="none"),
               design=list(galpha=1,gbeta=1,s=1,p=0),
               stvalue=list(s=1.8),upper=list(s=3),lower=list(s=1))

#Paphies donacina model in Table 4 of Cranfield et al (1996) with
#asymptotic model
data(P.donacina)
grotagplus(P.donacina,alpha=50,beta=80,
       model=list(mean="asymptotic",var="linear",seas="none"),
       design=list(galpha=1,gbeta=1,Lstar=0,s=1,nu=0,m=0,p=0),
       stvalue=list(galpha=10,gbeta=1.5,s=2),
       upper=list(galpha=15,gbeta=2.7,s=4),
       lower=list(galpha=7,gbeta=0.2,s=0.5),
       fixvalue=list(Lstar=80))

#Paphies donacina model in Table 4 of Cranfield et al (1996) with
#asymptotic model
data(P.donacina)
grotagplus(P.donacina,alpha=50,beta=80,
       model=list(mean="asymptotic",var="linear",seas="none"),
       design=list(galpha=1,gbeta=1,Lstar=0,s=1,nu=0,m=0,p=0),
       stvalue=list(galpha=10,gbeta=1.5,s=2),
       upper=list(galpha=15,gbeta=2.7,s=4),
       lower=list(galpha=7,gbeta=0.2,s=0.5),
       fixvalue=list(Lstar=80))

# Model 4 fit from Francis and Francis (1992) with different growth by sex
data(rig)
grotagplus(rig,dataID="Sex",alpha=70,beta=100,
           model=list(mean="Francis",var="linear",seas="none"),
          design=list(galpha=list("F","M"),gbeta=list("F","M"),s=1,nu=1,m=0,p=0),
          stvalue=list(galpha=c(5,4),gbeta=c(3,2),s=2,nu=0.5),
          upper=list(galpha=c(8,6),gbeta=c(5,4),s=4,nu=1),
          lower=list(galpha=c(3,2),gbeta=c(1.5,1),s=0.5,nu=0.2))

#Example where all parameters are fixed 
# to the values estimated values for model 4 of Francis and Francis (1992)]
grotagplus(rig,dataID="Sex",alpha=70,beta=100,
          model=list(mean="Francis",var="linear",seas="none"),
          design=list(galpha=0,gbeta=0,s=0,nu=0,m=0,p=0),
          stvalue=list(),upper=list(),lower=list(),
          fixvalue=list(galpha=list(design=list("F","M"),value=c(5.87,3.67)),
          gbeta=list(design=list("F","M"),value=c(2.52,1.73)),s=1.57,nu=0.58))
}
\keyword{misc}