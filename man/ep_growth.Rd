\name{ep_growth}
\alias{ep_growth}
\title{Fitting the von Bertalanffy growth model to length-stratified age samples}
\description{
  Estimation of von Bertanffy growth parameters based on length-stratified age samples (Perrault et al., 2020)
}
\usage{
ep_growth(len=NULL,age=NULL,Nh=NULL,nh=NULL,starts=list(Linf=60,
k=0.1,a0=-0.01,CV=0.5), 
bin_size=2,nlminb.control=list(eval.max=5000, 
iter.max=5000,trace=10),
tmb.control=list(maxit=5000,trace=FALSE),plot=TRUE)
}
\arguments{
  \item{len}{vector of lengths.}
  \item{age}{the vector of ages associated with the length vector.}
  \item{Nh}{the total sample size per bin. Includes the unaged fish.}
  \item{nh}{the total aged sample size per bin.}
  \item{starts}{the starting values for \emph{L-infinity}, \emph{K}, \emph{a0} and \emph{CV}. Required.}
  \item{bin_size}{the bin size (e.g., 2 for 2-cm) of the length stratification.}
  \item{nlminb.control}{controls for the \emph{nlminb} function. See function \emph{nlminb} for more information.}
  \item{tmb.control}{controls for the \emph{TMB} function. See package \emph{TMB} for more information.}
  \item{plot}{plot observed and model predicted lengths at age. Default is TRUE.}
}
\details{
 The von Bertalanffy growth model \emph{Lage=Linf*(1-exp(-K*(age-a0))} is fitted to length-at-age data 
 collected via length-stratified sampling following the EP method of Perreault et al. (2020).
 A plot of  model fit is generated unless \emph{plot}=FALSE. 
   
}
\value{
  List containing list elements of the model convergence, parameter estimates and predicted values.
}

\references{
Perrault, A. M. J., N. Zhang and Noel G. Cadigan. 2020. Estimation of growth parameters based on 
length-stratified age samples. Canadian Journal of Fisheries and Aquatic Sciences 77: 439-450. 
}
\author{Andrea Perrault,
Marine Institute of Memorial University of Newfoundland

 \email{andrea.perrault@mi.mun.ca}
 
}
\examples{
 \dontrun{
	data(ep.data)
	ep_growth(len=ep.data$length, age=ep.data$age, Nh=ep.data$Nh, nh=ep.data$nh,
        starts=list(Linf=60,
        k=0.1, CV=0.5 ,a0=-0.01), bin_size=2,
          nlminb.control=list(eval.max=5000, iter.max=5000, trace=10), 
          tmb.control=list(maxit=5000, trace=F),plot=TRUE)
 }
}
\keyword{misc}

