\name{growthlrt}
\alias{growthlrt}
\title{Likelihood Ratio Tests for Comparing Multiple Growth Curves}
\description{
  Likelihood ratio tests for comparison of two or more growth curves (von Bertalanffy, Gompertz and logistic)
}
\usage{
growthlrt(len = NULL, age = NULL, group = NULL, model = 1, error = 1,
 select = 1, Linf = c(NULL), K = c(NULL), t0 = c(NULL),plottype=0,
control=list(maxiter=10000,minFactor=1/1024,tol=1e-5))
}

\arguments{
  \item{len}{the vector of lengths of individual fish.}
  \item{age}{the vector of ages associated with the length vector.}
  \item{group}{the vector of character names specifying group association. The first character in the name must be a letter.}
  \item{model}{code indicating the growth model to use. 1 = von Bertalanffy, 2= Gompertz and 3= logistic. Default=1.}

  \item{error}{the error variance assumption.  1= constant variance for all \emph{lij}s; 2= constant variance for all mean 
lengths at age; 3=var of \emph{lij} varies with age. See methods a-c in Kimura (1980: pp. 766).  The required statistics 
for each type of error are calculated from the individual length-age observations.}
  \item{select}{the selection of starting values of \emph{L-infinity}, \emph{K}, and \emph{t0}. 1=automatic selection, 
2=user-specified. If \emph{select}=1, initial starting values of \emph{L-infinity}, \emph{K}, and \emph{t0} are 
 calculated from Walford lines (Everhart et al. 1975), and ages represented as decimal values are truncated to the 
integer before linear regression is applied. If select=2, the user must specify the values of \emph{L-infinity}, 
\emph{K}, and \emph{t0}.}
 \item{Linf}{if \emph{select}=2, the starting values of \emph{L-infinity} of the von Bertalanffy equation for each group.}
 \item{K}{if \emph{select}=2, the starting values of \emph{K} of the von Bertalanffy equation for each group.}
 \item{t0}{if \emph{select}=2, the starting values of \emph{t0} of the von Bertalanffy equation for each group.}
 \item{plottype}{the type of plot for each model. 1= observed versus predicted, 2= residuals. Default= 0 (no plot).}
 \item{control}{see function \emph{nls}.}
}
\details{
  Following Kimura (1980), the general model (one \emph{L-infinity}, \emph{K}, and \emph{t0} for each group) 
  and four sub models are fitted to the length and age data using function \emph{nls} (nonlinear least squares). For
 each general model-sub model comparison,   likelihood ratios are calculated by using the residual sum-of-squares and 
are tested against chi-square statistics with the appropriate degrees of freedom. Individual observations of 
lengths-at-age are required. If error variance assumptions 2 or 3, mean lengths and required statistics are calculated. 
The parameters are fitted using a model.matrix where the 1st column is a row of 1s representing
the parameter estimate of the reference group (lowest alpha-numeric order) and the remaining group columns
have 1 if group identifier is the current group and 0 otherwise. The group number depends on the alph-numeric order.
See function \emph{model.matrix}.  

The model choices are: 

von Bertalanffy
La=Linf(1-exp(-K*(a-t0)))

Gompertz
La=Linf*exp(-exp(-K*(a-t0)))

Logisitic
La=Linf/(1+exp(-K*(a-t0)))
 

To extract the growth parameters for each group under an hypothesis: 
   
            x$'model Ho'$coefficients

            x$'model H1'$coefficients

            x$'model H2'$coefficients

            x$'model H3'$coefficients

            x$'model H4'$coefficients

where \emph{x} is the output object. 

As an example, let's say three groups were compared.To get the L-infinity estimates for each groups,

Linf1<-x$'model Ho'$coefficients[1]

Linf2<-Linf1+ x$'model Ho'$coefficients[2]

Linf3<-Linf1+ x$'model Ho'$coefficients[3]

For models H1, H2, H3 and H4, the parameter L1 or K1 or t01 will be
shared across groups. 

If RSSHX >RSSH0, less information is accounted for by RSSHX model (where X is hypothesis 1, 2,..etc.).  
If Chi-square is significant, RSSH0 is the better model. If Chi-square is not significant, RSSHX is 
the better model. 
}
\value{
  \item{results }{list element with the likelihood ratio tests comparing von Bertalanffy models.}
  \item{model Ho}{list element with the \code{nls} fit for the general model.}
  \item{model H1}{list element with the \code{nls}  for model H1 (Linf1=Linf2=..=Linfn) where n is the number of groups.}
  \item{model H2}{list element with the \code{nls} fit for model H2 (K1=K2=..=Kn).}
  \item{model H3}{list element with the \code{nls} fit for model H3 (t01=t02=...=t0n).}
  \item{model H4}{list element with the \code{nls} fit for model H4 (Linf1=Linf2=..=Linfn, K1=K2=..=Kn, t01=t02=...=t0n).}
  \item{rss}{list element with the residual sum-of-squares from each model.}
  \item{residuals}{list element with the residuals from each model.}
}
\references{Everhart, W. H., A. W. Eipper, and W. D. Youngs. 1975. Principles of Fishery Science. Cornell 
University Press.

Kimura, D. K. 1980. Likelihood methods for the von Bertalanffy growth curve. U. S. Fish. Bull. 77(4): 765-776.  
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
## Normally, the length and age data will represent data for individuals.  
## Kimura's data are mean lengths-at-age but are usable because error=2 
## will calculate mean lengths-at-age from individual data. Since only  
## one value is present for each age,the mean length will be calculated
## as the same value.
data(Kimura)
growthlrt(len=Kimura$length,age=Kimura$age,group=Kimura$sex,model=1,error=2,select=1,
plottype=2)
}
\keyword{misc}
\keyword{htest}
