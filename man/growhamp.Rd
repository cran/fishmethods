\name{growhamp}
\alias{growhamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
von Bertalanffy Growth Models for Tagging Data Incorporating Individual Variation  
}
\description{
Function fits growth models of Hampton (1991) to length and time-at-large data from tagging studies 
}
\usage{
growhamp(L1 = NULL, L2 = NULL, TAL = NULL, 
models = c(1, 2, 3, 4, 5, 6, 7), 
method = c("Nelder-Mead", "Nelder-Mead", "Nelder-Mead",
 "Nelder-Mead", "Nelder-Mead", "Nelder-Mead", "Nelder-Mead"), 
varcov = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), 
Linf = list(startLinf = NULL, lowerLinf = NULL, upperLinf = NULL), 
K = list(startK = NULL, lowerK = NULL, upperK = NULL), 
sigma2_error = list(startsigma2 = NULL, lowersigma2 = NULL, uppersigma2 = NULL), 
sigma2_Linf = list(startsigma2 = NULL, lowersigma2 = NULL, uppersigma2 = NULL),
sigma2_K = list(startsigma2 = NULL, lowersigma2 = NULL, uppersigma2 = NULL), 
mu_measure = 0, sigma2_measure = 0, 
control = list(maxit = 1000))
}
\arguments{
  \item{L1}{
  Vector of release lengths. Each row presents the length of an individual.
}
  \item{L2}{
Vector of recapture lengths.
}
  \item{TAL}{
  vector of associated time-at-large data. Calculated as the recapture date minus release date.
}
  \item{models}{
 The models to fit. 1 = Faber model, 2 = Kirkwood and Somers model, 3 = Kirkwood and Somers model with model error, 4 = Kirkwood and Somers model with model and release-length-measurement error, 5 = Sainsbury model, 6 = Sainsbury model with model error, and 7 = Sainsbury model with model and release-length-measurement error.  Default is all: c(1,2,3,4,5,6,7).
}
  \item{method}{
Character vector of optimization methods used in \code{optim} to solve parameters for each model. A different method can be selected for each model. Choices are "Nelder-Mead","BFGS","CG","L-BFGS-B" and "SANN". See help for \code{optim}. Default is "Nelder-Mead". If there are fewer values specified in \code{method} than the number specified in \code{models}, a warning message is produced and the last value in the \code{method} vector is used for the remaining models.  
}
  \item{varcov}{
  Logical vector specifying whether the parameter variance-covariance matrix of each model should be outputted. A different logical can specified for each model. If there are fewer values specified in \code{varcov} than the number specified in \code{models}, a warning message is produced and the last value in the \code{varcov} vector is used for the remaining models.  
}
  \item{Linf}{
A list of starting (startLinf), lower bound (lowerLinf) and upper bound (upperLinf) of Linfinity of the von Bertalanffy equation used in the optimization. The lower and upper bounds are used only with method "L-BFGS-B".
}
  \item{K}{
A list of starting (startK), lower bound (lowerK) and upper bound (upperK) of K (growth coefficient) of the von Bertalanffy equation used in the optimization.
The lower and upper bounds are used only with method "L-BFGS-B".
}
  \item{sigma2_error}{
A list of starting (startsigma2), lower bound (lowersigma2) and upper bound (uppersigma2) of the error variance used in the optimization. The lower and upper bounds are used only with method "L-BFGS-B". This parameter is used in models 1,3,4,6 and 7.
}
  \item{sigma2_Linf}{
A list of starting (startsigma2), lower bound (lowersigma2) and upper bound (uppersigma2) of the Linfinity variance used in the optimization. The lower and upper bounds are used only with method "L-BFGS-B". This parameter is used in models 2,3,4,5,6,and 7.
}
  \item{sigma2_K}{
A list of starting (startsigma2), lower bound (lowersigma2) and upper bound (uppersigma2) of the K (growth coefficient) variance used in the optimization. The lower and upper bounds are used only with method "L-BFGS-B". This parameter is used in models 5,6, and 7.
}
  \item{mu_measure}{
Release measurement error. This parameter is used in models 4 and 7.
Default=0.
}
  \item{sigma2_measure}{
Variance of release measurement error.  This parameter is used in models 4 and 7.
Default=0.
}
  \item{control}{
A list of control parameters for \code{optim}. See function \code{optim} for details.
}
}
\details{
The seven models are fitted by maximum likelihood using formulae shown in Hampton 1991. Due to the number of parameters estimated, some models can be sensitive to the initial starting values. It is recommended that the starting values are tested for sensitivity to ensure the global minimum has been reached. Sometimes, the hessian matrix, which is inverted to obtain the variance-covariance matrix, will not be positive, definite and therefore will produce an error. Again, try different starting values for parameters and lower and upper bounds if applicable. 
}
\value{
 \item{results}{list element containing the parameter estimates in table format for each model. Column names are \code{model}, \code{Linf}, \code{K}, \code{s2Linf} (variance of Linf), \code{s2K} (variance of K), \code{s2error} (error variance), \code{boundary} (0 = no issues; 1 = one or more parameter estimates are at constraint boundaries), \code{-Log Likelihood}, \code{AIC} (Akaike's Information Criterion, and \code{method}
}
\item{varcov}{if varcov=TRUE, list element containing the variance-covariance matrix for each model.}  
\item{residuals}{list element containing the residuals (observed-predicted values) for each model.}
}
\references{Hampton, J. 1991. Estimation of southern bluefin tuna Thunnus maccoyii growth parameters from tagging data, using von Bertalanffy models incorporating individual variation. U. S. Fishery Bulletin 89: 577-590.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{mort.al}}
}
\examples{
\dontrun{
## Models 1,2 and 3 below are models 1,2, and 4 in Table 4.17 of ##Quinn and Deriso 
data(trout)
growhamp(L1=trout$L1,L2=trout$L2,TAL=trout$dt,models=c(1,2,3),
       method=c("Nelder-Mead","Nelder-Mead","L-BFGS-B"),
       varcov=c(TRUE,TRUE,TRUE),
       Linf=list(startLinf=650,lowerLinf=400,upperLinf=800),       
       K=list(startK=0.30,lowerK=0.01,upperK=1),
       sigma2_error=list(startsigma2=100,lowersigma2=0.1,uppersigma2=10000),
       sigma2_Linf=list(startsigma2=100,lowersigma2=0.1,uppersigma2=100000),	
       sigma2_K=list(startsigma2=0.5,lowersigma2=1e-8,uppersigma2=10))
}
}
\keyword{misc}
