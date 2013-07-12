\name{grotag}
\alias{grotag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum likelihood estimation of growth and growth variability from tagging data - Francis (1988)}
\description{
This function estimates parameters of Francis (1988)'s growth model
using tagging data. The data are fitted using a constrained maximum
likelihood optimization performed by optim using the "L-BFGS-B" method.
}
\usage{
grotag(L1 = NULL, L2 = NULL, T1 = NULL, T2 = NULL, alpha = NULL, beta = NULL, 
 design = list(nu = 0, m = 0, p = 0, sea = 0), 
 stvalue = list(sigma = 0.9, nu = 0.4, m = -1, p = 0.01, u = 0.4, w = 0.4), 
 upper = list(sigma = 5, nu = 1, m = 2, p = 1, u = 1, w = 1), 
 lower = list(sigma = 0, nu = 0, m = -2, p = 0, u = 0, w = 0), gestimate = TRUE, 
 st.ga = NULL, st.gb = NULL, st.galow = NULL, st.gaup = NULL, st.gblow = NULL,
 st.gbup = NULL, control = list(maxit = 10000))
}
\arguments{
  \item{L1}{Vector of length at release of tagged fish}
  \item{L2}{Vector of length at recovery of tagged fish}
  \item{T1}{Vector of time at release of tagged fish}
  \item{T2}{Vector of time at recovery of tagged fish}
  \item{alpha}{Numeric value giving an arbitrary length alpha}
  \item{beta}{Numeric value giving an arbitrary length beta (\code{beta} > \code{alpha})}
  \item{design}{List specifying the design of the model to estimate. Use 1 to designate whether 
       a parameter(s) should be estimated. Type of parameters are: nu=growth variability (1 parameter), 
       m=bias parameter of measurement error (1 parameter), p=outlier probability (1 parameter), 
       and sea=seasonal variation (2 parameters: u and w).  Model 1 of Francis is the default 
       settings of 0 for nu, m, p and sea.}
  \item{stvalue}{Starting values of sigma (s) and depending on
    the design argument, nu, m, p, u, and w used as input in the nonlinear estimation
       (function \emph{optim}) routine.}
  \item{upper}{Upper limit of the model parameters' (nu, m, p, u, and w) region to be investigated.}
  \item{lower}{Lower limit of the model parameters' (nu, m, p, u, and w) region to be investigated.}
  \item{gestimate}{Logical specifying whether starting values of ga and gb (growth increments of alpha and beta) 
      should be estimated automatically.  Default = TRUE.}
  \item{st.ga}{If gestimate=FALSE, user-specified starting value for ga.}
  \item{st.gb}{If gestimate=FALSE, user-specified starting value for gb.}
  \item{st.galow}{If gestimate=FALSE, user-specified lower limit for st.ga used in optimization.}
  \item{st.gaup}{If gestimate=FALSE, user-specified upper limit for st.ga used in optimization.}
  \item{st.gblow}{If gestimate=FALSE, user-specified lower limit for st.gb used in optimization.}
  \item{st.gbup}{If gestimate=FALSE, user-specified upper limit for st.gb used in optimization.}
   \item{control}{Additional controls passed to the optimization function \emph{optim}.}
}
\details{
The methods of Francis (1988) are used on tagging data to the estimate of growth and growth variability.  
The estimation of all models discussed is allowed.  The growth variability defined by equation 5 
in the reference is used throughout.  
}
\value{
  \item{table}{list element containing the model output similar to Table 3 of Francis (1988).  
         The Akaike's Information Criterion (AIC) is also added to the output.}
  \item{VBparms}{list element containing the conventional paramaters of the von
         Bertalanffy model (Linf and K).}  
  \item{correlation}{list element containing the parameter correlation matrix.}
  \item{predicted}{list element containing the predicted values from the model.}
  \item{residuals}{list element containing the residuals of the model fit.}
}
\references{Francis, R.I.C.C., 1988. Maximum likelihood estimation of growth and growth variability 
from tagging data. New Zealand Journal of Marine and Freshwater Research, 22, p.42-51. 
}
\author{
Marco Kienzle \email{Marco.Kienzle@gmail.com}

Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}

\examples{
data(bonito)

#Model 4 of Francis (1988)
with(bonito,
 grotag(L1=L1, L2=L2, T1=T1, T2=T2,alpha=35,beta=55,
 	design=list(nu=1,m=1,p=1,sea=1),
 	stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0,u=0.4,w=0.4),
 	upper=list(sigma=5,nu=1,m=2,p=0.5,u=1,w=1),
 	lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=1e4)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
