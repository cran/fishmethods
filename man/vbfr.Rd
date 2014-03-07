\name{vbfr}
\alias{vbfr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Francis' re-parameterization of the von Bertalanffy growth equation for length-age data}
\description{
Fits the re-parameterized von Bertalanffy growth equation of Francis (1988) 
by using nonlinear least-squares}
\usage{
vbfr(age = NULL, L = NULL, agephi = NULL, agepsi = NULL, graph = TRUE, 
gestimate = TRUE, Lphiparms = c(NA, NA, NA), Lchiparms = c(NA, NA, NA), 
Lpsiparms = c(NA, NA, NA),control = list(maxiter = 10000))
}
\arguments{
  \item{age}{Vector of ages of individual fish.}
  \item{L}{Vector of lengths of individual fish.}
  \item{agephi}{Arbitrary reference age phi}
  \item{agepsi}{Arbitrary reference age psi. agepsi>agephi.}
  \item{graph}{Logical specifiying whether observed versus predicted, and 
    residual plots should be drawn. Default=TRUE.}
  \item{gestimate}{Logical specifying whether automatic generation of starting
       values of \emph{lphi}, \emph{lchi} and \emph{lpsi} should be used. Default=TRUE. 
       If gestimate=FALSE, user-specified starting, lower and upper limits of parameters
       must be entered.}
  \item{Lphiparms}{If gestimate=FALSE, starting value, lower limit and upper limit 
       of \emph{lphi} used in \emph{nls}.}
  \item{Lchiparms}{If gestimate=FALSE, starting value, lower limit and upper limit 
       of \emph{lchi} used in \emph{nls}.}
  \item{Lpsiparms}{If gestimate=FALSE, starting value, lower limit and upper limit 
       of \emph{lpsi} used in \emph{nls}.}
  \item{control}{see \code{control} under function \emph{nls}.}
}
\details{
Francis (1988) re-parameterized the von Bertalanffy growth equation for age-length 
in order to make equivalent comparison of parameters to parameters of a common model 
used to estimate growth from tagging data.  Three parameters,  \emph{lphi}, \emph{lchi} 
and \emph{lpsi}, are estimated. The re-parameterization also has better 
statistical properties than the original equation.

The formulae to get the conventional von Bertalanffy parameters are: 
 
Linf = lphi + (lpsi-lphi)/(1-r^2) where r = (lpsi-lchi)/(lchi-lphi)

K = -(2*log(r))/(agepsi-agephi)

t0 = agephi + (1/K)*log((Linf-lphi)/Linf)

If gestimate=TRUE, unconstrained nonlinear least-squares (function \emph{nls}) is used 
   to fit the model.  If gestimate=FALSE, constrained nonlinear least-squares 
   is used (algorithm "port" in \emph{nls}). 
}
\value{
\emph{nls} object of model results. Use \emph{summary} to extract results.
}
\references{
Francis, R. I. C. C. 1988. Are growth parameters estimated from tagging and age-length data comparable?
Can. J. Fish. Aquat. Sci. 45: 936-942.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}
}
\examples{
data(pinfish)
with(pinfish,vbfr(age=age,L=sl,agephi=3,agepsi=6))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
