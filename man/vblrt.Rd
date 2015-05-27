\name{vblrt}
\alias{vblrt}
\title{Likelihood Ratio Tests for Comparing Two von Bertalanffy Growth Curves}
\description{
  Likelihood ratio tests for comparison of two von Bertalanffy growth curves following Kimura (1980).
}
\usage{
vblrt(len = NULL, age = NULL, group = NULL, error = 1,
 select = 1, Linf = NULL, K = NULL, t0 = NULL,plottype=0,
control=list(maxiter=10000,minFactor=1/1024,tol=1e-5))
}
\arguments{
  \item{len}{the vector of lengths of individual fish.}
  \item{age}{the vector of ages associated with the length vector.}
  \item{group}{the vector of character or numeric codes specifying group association.}
  \item{error}{the error variance assumption.  1= constant variance for all \emph{lij}s; 2= constant variance for all mean 
lengths at age; 3=var of \emph{lij} varies with age. See methods a-c in Kimura (1980: pp. 766).  The required statistics 
for each type of error are calculated from the individual length-age observations.}
  \item{select}{the selection of starting values of \emph{L-infinity}, \emph{K}, and \emph{t0}. 1=automatic selection, 
2=user-specified. If \emph{select}=1, initial starting values of \emph{L-infinity}, \emph{K}, and \emph{t0} are 
 calculated from Walford lines (Everhart et al. 1975), and ages represented as decimal values are truncated to the 
integer before linear regression is applied. If select=2, the user must specify the values of \emph{L-infinity}, 
\emph{K}, and \emph{t0}.}
 \item{Linf}{if \emph{select}=2, the starting value for \emph{L-infinity} of the von Bertalanffy equation.}
 \item{K}{if \emph{select}=2, the starting value for \emph{K} of the von Bertalanffy equation.}
 \item{t0}{if \emph{select}=2, the starting value for \emph{t0} of the von Bertalanffy equation.}
 \item{plottype}{the type of plot for each model. 1= observed versus predicted, 2= residuals. Default= 0 (no plot).}
 \item{control}{see function \emph{nls}.}
}
\details{
  Following Kimura (1980), the general model (6 parameters; one \emph{L-infinity}, \emph{K}, and \emph{t0} for each group) 
  and four sub models are   fitted to the length and age data using function \emph{nls} (nonlinear least squares). For
 each general model-sub model comparison,   likelihood ratios are calculated by using the residual sum-of-squares and 
are tested against chi-square statistics with the appropriate degrees of freedom. Individual observations of 
lengths-at-age are required. If error variance assumptions 2 or 3, mean lengths and required statistics are calculated. 
A dummy vector called \emph{cat}, containing 0s for the first group with lower alpha-numeric order and 1s 
  for the second group, is used in the estimation of group parameters. To extract the VB parameters for each group 
  from the Ho model:
            # Group 1
            Linf1<-x$'model Ho'$coefficients[1]
		K1<-x$'model Ho'$coefficients[3]
		t01<-x$'model Ho'$coefficients[5]

            # Group 2
		Linf2<-Linf1+x$'model Ho'$coefficients[2]
		K2<-K1+x$'model Ho'$coefficients[4]
		t02<-t01+x$'model Ho'$coefficients[6]
	
where \emph{x} is the output object. For models H1, H2, H3 and H4, the shared VB parameter for the second group 
will be the same as the first group. 
}
\value{
  \item{results }{list element with the likelihood ratio tests comparing von Bertalanffy models.}
  \item{model Ho}{list element with the \code{nls} fit for the general model.}
  \item{model H1}{list element with the \code{nls}  for model H1 (Linf1=Linf2).}
  \item{model H2}{list element with the \code{nls} fit for model H2 (K1=K2).}
  \item{model H3}{list element with the \code{nls} fit for model H3 (t01=t02).}
  \item{model H4}{list element with the \code{nls} fit for model H4 (Linf1=Linf2, K1=K2, t01=t02).}
  \item{rss}{list element with the residual sum-of-squares from each model.}
  \item{residuals}{list element with the residuals from each model.}
}
\references{Everhart, W. H., A. W. Eipper, and W. D. Youngs. 1975. Principles of Fishery Science. Cornell 
University Press.

Kimura, D. K. 1980. Likelihood methods for the von Bertalanffy growth curve. U. S. Fish. Bull. 77(4): 765-776.  
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
## Normally, the length and age data will represent data for individuals.  
## Kimura's data are mean lengths-at-age but are usable because error=2 
## will calculate mean lengths-at-age from individual data. Since only  
## one value is present for each age,the mean length will be calculated
## as the same value.
data(Kimura)
vblrt(len=Kimura$length,age=Kimura$age,group=Kimura$sex,error=2,select=1,
plottype=2)
}
\keyword{misc}
\keyword{htest}
