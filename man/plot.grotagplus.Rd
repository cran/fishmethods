\name{plot.grotagplus}
\alias{plot.grotagplus}
\title{
Plotting Tagging-Growth Objects
}
\description{
 Plotting method for output from function grotagplus, which has
  class "grotagplus".
}
\usage{
\method{plot}{grotagplus}(x,plot.type="meangrowth",Linitial=NULL,resid.spec=list(Pearson=T,
      x="mean.delL"),xlim=NULL,ylim=NULL,pch=20,leg.loc=NULL,
      age.based.growth=NULL,\dots)
}
\arguments{
    \item{x}{Growth-model fit to tagging data as output by function
    "grotagplus".}
  \item{plot.type}{Character string identifying the type of plot
    required: "meangrowth" = mean annual growth vs initial length;
    "traj" = one-year growth trajectory of fish of initial
    length specified by \code{Linitial}; or "resid" = plot of ordinary
    or Pearson residuals (plot details specified by \code{resid.spec}).}
  \item{Linitial}{Initial length to use for plot of growth trajectory.}
  \item{resid.spec}{List, specifying details of a residual plot, with
    components "Pearson" (logical, if T [default] plot Pearson residuals, 
    otherwise simple residuals) and "x" (the x-variable in the plot - either
    "L1", length at tagging; "delT", time at liberty; or
    "mean.delL", expected length increment).}
  \item{xlim}{Allow the user to set x-limits for a plot
    that differ from those defined by the range of the plotted data.}
   \item{ylim}{Allow the user to set y-limits for a plot
    that differ from those defined by the range of the plotted data.}
  \item{pch}{Allows the user to change the plotting symbol for residual
    plots from the default pch=20.}
  \item{leg.loc}{Allows the user to change the legend location
    from its default position ("topright" for meangrowth and resid;
    "topleft" for traj).  Note that a legend is used only for traj or
    for other plots with multiple datasets.}
  \item{age.based.growth}{This argument allows the user to add,
    to a meangrowth plot, growth estimates (plotted as dashed lines)
    from age-length datasets.  It should be a list of vectors,
    each of which contains estimates of mean length corresponding to a
    vector of increasing ages whose increments are  always 1 year
    (the ages are not included in the argument because they are not used
    in the plot, and the age vectors need not be the same in each
    component).  If the list is named then the names will be interpreted
    as identifying different datasets.  If a name appears in
    fit$datasetnames the age-based growth will be plotted with the same
    colour as the corresponding tagging growth.  If the list is not
    named then it must be of the same length as fit$datasetnames
    (or of length 1 if there is only one dataset in the tagging data)
    and it will be assumed that the ith component corresponds to the
    ith tagging dataset.} 
\item{...}{Other graphical parameters. See \code{par}}
}
\details{
  Examples of the three plot types are given in  Figs 7 & 8 of
  Francis and Francis (1992), for "resid" and "meangrowth",
  respectively; and in Fig. 2 of Francis (1988), for "traj". 

  plot.type="meangrowth" is the recommended way for plotting growth
  rates estimated from tagging data.  Argument age.based.growth allows a
  rough comparison between these growth estimates and those from
  age-length data (the comparison is between the mean growth at length L
  and that at the age for which the mean length is L).
    
  The traj plot, as well as showing the mean (i.e., expected) growth (solid
  line), shows 95% confidence intervals for growth, both without
  (dashed lines) and with (dotted lines) allowance for measurement
  error.

  In residual plots, a dashed lowess line is plotted for each dataset
  to indicate any trend and, for Pearson residuals, dotted lines at +/- 2 
  indicate approximate 95% confidence intervals.
  
  For fits using multiple datasets, colour is used to distinguish the
  datasets.  Use "palette" to change the match between colour and
  dataset (the ith colour in the palette is associated with the ith
  element in fit$datasetnames).
}

\references{
\code{1} Francis, R.I.C.C., 1988. Maximum likelihood estimation of
  growth and growth variability from tagging data.
  New Zealand Journal of Marine and Freshwater Research, 22, p.42-51.
  
\code{2} Francis, M.P. and Francis, R.I.C.C. 1992.  Growth rate
  estimates for New Zealand rig (Mustelus lenticulatus).  Australian
  Journal of Marine and Freshwater Research 43: 1157-1176
}
\author{
Chris Francis \email{chrisfrancis341@gmail.com}

Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}

Marco Kienzle \email{Marco.Kienzle@gmail.com}
}
\seealso{\code{\link{grotagplus}} \code{\link{print.grotagplus}}}

\examples{
# Plot of mean growth like that in Fig 8. of Francis & Francis (1992)
data(rig)
fit <- grotagplus(rig,dataID="Sex",alpha=70,beta=100,
                  model=list(mean="Francis",var="linear",seas="none"),
                 design=list(galpha=list("F","M"),gbeta=list("F","M"),
                             s=1,nu=1,m=0,p=0),
                  stvalue=list(galpha=c(5,4),gbeta=c(3,2),s=2,nu=0.5),
                  upper=list(galpha=c(8,6),gbeta=c(5,4),s=4,nu=1),
                 lower=list(galpha=c(3,2),gbeta=c(1.5,1),s=0.5,nu=0.2))
mnlenatage <- list(F=90.7*(1-exp(-0.42*(seq(1.5,6.5)-0.77))),
           M= 118.7*(1-exp(-0.16*(seq(4,11)-2.02))),
           PGM=161.1*(1-exp(-0.11*(seq(3.5,10.5)-1.91))))
plot(fit,age.based.growth=mnlenatage)
## Residual plots
fit <- grotagplus(rig,dataID="Sex",alpha=70,beta=100,
                  model=list(mean="Francis",var="linear",seas="none"),
                 design=list(galpha=list("F","M"),gbeta=list("F","M"),
                             s=1,nu=1,m=0,p=0),
                  stvalue=list(galpha=c(5,4),gbeta=c(3,2),s=2,nu=0.5),
                  upper=list(galpha=c(8,6),gbeta=c(5,4),s=4,nu=1),
                 lower=list(galpha=c(3,2),gbeta=c(1.5,1),s=0.5,nu=0.2))
plot(fit,"resid")
plot(fit,"resid",resid.spec=list(Pearson=FALSE,x="L1"))
## Trajectory plot as in Fig. 2 of Francis (1988)
data(bonito)
fit <- grotagplus(bonito,alpha=35,beta=55,
               design=list(galpha=1,gbeta=1,s=1,nu=1,m=1,p=1,u=1,w=1),
               stvalue=list(s=0.81,nu=0.3,m=0,p=0.01,u=0.5,w=0.5),
               upper=list(s=3,nu=1,m=2,p=0.1,u=1,w=1),
               lower=list(s=0.1,nu=0.1,m=-2,p=0,u=0,w=0))
plot(fit,"traj",Linitial=35)
}
\keyword{misc}