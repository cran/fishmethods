\name{print.grotagplus}
\alias{print.grotagplus}
\title{
Printing Tagging-Growth Objects
}
\description{
 Printing method for output from function grotagplus, which has
  class "grotagplus". 
}
\usage{
\method{print}{grotagplus}(x,precision=c(est="sig3",stats="dec1",cor="dec2"),\dots)
}
\arguments{
   \item{x}{Growth-model fit to tagging data as output by function
    "grotagplus".}
  \item{precision}{Named character vector specifying the printing
    precision for each of three categories of output: "est" (applies
    to fixed and estimated parameters and to Linf.k); "stats" (for
    negloglikl and AIC); and "cor" (for the parameter correlation
    matrix).  Values should be either "sigx", for x significant figures, or
    "decx" for x decimal places.}
  \item{...}{Other print parameters.}
}
\details{
 Outputs from grotagplus are produced to a precision which is
  usually much greater than is warranted.  To see this full precision
  print individual components, e.g., print(fit$parest).
}

\author{
Chris Francis \email{chrisfrancis341@gmail.com}

Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}

Marco Kienzle \email{Marco.Kienzle@gmail.com}
}
\seealso{\code{\link{grotagplus}} \code{\link{plot.grotagplus}}}

\examples{
#Model 4 of Francis (1988)
data(bonito)
fit <- grotagplus(bonito,alpha=35,beta=55,
               design=list(galpha=1,gbeta=1,s=1,nu=1,m=1,p=1,u=1,w=1),
               stvalue=list(s=0.81,nu=0.3,m=0,p=0.01,u=0.5,w=0.5),
               upper=list(s=3,nu=1,m=2,p=0.1,u=1,w=1),
               lower=list(s=0.1,nu=0.1,m=-2,p=0,u=0,w=0))
print(fit)
}
\keyword{misc}