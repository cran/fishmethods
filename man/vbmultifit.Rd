\name{vbmultifit}
\alias{vbmultifit}
\title{Fit a Multi-Group von Bertalanffy Growth Model}
\description{
  Fits a von Bertalanffy growth curve to length and age for two or more groups.
}
\usage{
vbmultifit(len=NULL,age=NULL,group=NULL,fixed=c(1,1,1),error=1,
        select=1,Linf=c(NULL),K=c(NULL),t0=c(NULL),plot=FALSE,
                control=list(maxiter=10000,minFactor=1/1024,tol=1e-5))
}
\arguments{
  \item{len}{the vector of lengths of individual fish.}
  \item{age}{the vector of ages associated with the length vector.}
  \item{group}{the vector of character names specifying group association. The first character in the name must be a letter.}
  \item{fixed}{arguments specifying that Linf, K or t0 should be fitted as a constant 
  between groups or as separate parameters for each group. 1 = single parameter between 
  groups, 2 = separate parameters for each group. The order of \emph{fixed} is c(Linf,K,t0).} 
 \item{error}{the error variance assumption.  1= constant variance for all \emph{lij}s; 2= constant variance for all mean 
lengths at age; 3=var of \emph{lij} varies with age. See methods a-c in Kimura (1980: pp. 766).  The required statistics 
for each type of error are calculated from the individual length-age observations.}
  \item{select}{the selection of starting values of \emph{L-infinity}, \emph{K}, and \emph{t0}. 1=automatic selection, 
2=user-specified. If \emph{select}=1, initial starting values of \emph{L-infinity}, \emph{K}, and \emph{t0} are 
 calculated from Walford lines (Everhart et al. 1975), and ages represented as decimal values are truncated to the 
integer before linear regression is applied. If select=2, the user must specify values of \emph{L-infinity}, 
\emph{K}, and \emph{t0} for each group.}
 \item{Linf}{if \emph{select}=2, the starting values for \emph{L-infinity} of the von Bertalanffy equation, one for each group.}
 \item{K}{if \emph{select}=2, the starting values for \emph{K} of the von Bertalanffy equation, one for each group.}
 \item{t0}{if \emph{select}=2, the starting value for \emph{t0} of the von Bertalanffy equation, one for each group.}
 \item{plot}{logical argument specifying whether observed versus predicted and residuals graphs should be plotted. Default is FALSE.}
 \item{control}{see function \emph{nls}.}
}
\details{
 A single von Bertalanffy model is fitted to the length and age data of two or more groups using function \emph{nls} (nonlinear least squares). Parameters can be estimated for each group or as constants across groups. Individual observations of lengths-at-age are required. If error variance assumptions 2 or 3, mean lengths and required statistics are calculated. The parameters are fitted using a model.matrix where the 1st column is a row of 1s representing the parameter estimate of the reference group (group with lowest alpha-numeric order) and the remaining group columns have 1 if group identifier is the current group and 0 otherwise. See function \emph{model.matrix}.  
This is a companion function to function \emph{vblrt}.

When seperate parameters are estimated for each group, estimates for the the non-reference groups would be the reference-group estimated parameters (e.g., Linf1 or K1 or t01) plus the coefficent estimate for the nth group (e.g., group 2: Linf2 or K2, or t02) based on the alpha-numeric order.  If the parameter is assumed constant across groups, then estimates of Linf1 or K1 or t01 is used as the parameter for each group.
}
\value{
  \item{results}{list element containing summary statistics of \emph{nls} fit}
  \item{residuals}{list element with the residuals from the model.}
}
\references{Everhart, W. H., A. W. Eipper, and W. D. Youngs. 1975. Principles of Fishery Science. Cornell 
University Press.

Kimura, D. K. 1980. Likelihood methods for the von Bertalanffy growth curve. U. S. Fish. Bull. 77(4): 765-776.  
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{vblrt}}}
\examples{
data(Kimura)
vbmultifit(len=Kimura$length,age=Kimura$age,group=Kimura$sex,fixed=c(2,1,1),
error=1,select=1,Linf=NULL,K=NULL,t0=NULL,plot=FALSE,control=list(maxiter=10000,
minFactor=1/1024,tol=1e-5))
}
\keyword{misc}
