\name{growth}
\alias{growth}
\title{Fitting Growth Curves to Length- or Weight-at-Age Data}
\description{
  Fits various growth models to length and weight-at-age data.
}
\usage{
growth(intype=1,unit=1,size=NULL,age=NULL,calctype=1,wgtby=1,s2=NULL,error=1,
         Sinf=NULL,K=NULL,t0=NULL,B=3,graph=TRUE,
         control=list(maxiter=10000,minFactor=1/1024,tol=1e-5))
}
\arguments{
  \item{intype}{the input format: 1= individual size data; 2 = mean size data. Default intype=1.}
  \item{unit}{the size unit: 1= length; 2 = weight. Default unit=1.}
  \item{size}{the vector of size (length or weight) data.}
  \item{age}{the vector of ages associated with the size vector.}
  \item{calctype}{if intype=1, 1 = use individual size data; 2 = calculate mean size from individual size data. Default calctype=1.}
  \item{wgtby}{weighting scheme: 1 = no weighting; 2 = weight means by inverse variance. Weighting of individual data points is not allowed. Default wgtby=1.}
  \item{s2}{if intype=2 and wgtby=2, specify vector of sample variances associated with mean size-at-age data.}
  \item{error}{the error structure: 1 = additive; 2 = multiplicative. Default error=1.}
  \item{Sinf}{the starting value for \emph{L-infinity or W-infinity} of the growth models. Required.}
 \item{K}{the starting value for \emph{K} of the growth models.}
 \item{t0}{the starting value for \emph{t0} of the growth models.}
 \item{B}{the length-weight equation exponent used in the von Bertalanffy growth model for weight. Default B=3.}
 \item{graph}{logical value specifying if fit and residual plots should be drawn. Default graph = TRUE.}
 \item{control}{see function \emph{nls}.}
}
\details{
  If unit is length, a standard von Bertalanffy growth curve is fitted to length-at-age data using nonlinear least-squares (function \emph{nls}).
  If individual data are provided, mean size data can be calculated by specifying \emph{calctype}=2. When fitting mean size data, 
  observations can be weighted by the inverse sample variance(\emph{wgtby}=2), resulting in weighted nonlinear least squares. Additive or multiplicative 
  error structures are specified via \emph{error}.  See page 135 in Quinn and Deriso (1999) for more information on error structures.
 
  If unit is weight, three growth models (von Bertalanffy, Gompert and logistic) are fitted to weight-at-age data. 
  The exponent for the von Bertalanffy growth in weight model is not estimated and must be specified (\emph{B}).
  Individual and mean weight data can be used. When fitting mean size data, observations can be weighted by the inverse 
  sample variance (\emph{wgtby}=2), resulting in weighted nonlinear least squares. Additive or multiplicative 
  error structures are specified via \emph{error}.
 
  Plots of model fit and residuals are generated unless \emph{graph}=FALSE.   
}
\value{
  List containing list elements of the equation/structure and \emph{nls} output for each model.
  Information from \emph{nls} output can be extracted using standard functions (e.g., \emph{summary()}).
 }

\references{Quinn, T. J. and R. B. Deriso. 1999. Quantitative fish dynamics. Oxford University Press. 542 pages. 
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
data(pinfish)
growth(intype=1,unit=1,size=pinfish$sl,age=pinfish$age,
        calctype=1,wgtby=1,error=1,Sinf=200,K=0.3,t0=-1)
}
\keyword{misc}

