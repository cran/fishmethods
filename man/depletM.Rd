\name{depletM}
\alias{depletM}
\title{Catch-Effort Depletion Methods Incorporating Natural Mortality}
\description{
   Following Rosenberg et al. (1990), the Leslie-Davis and catch equation depletion methods that 
incorporate natural mortality are implemented for a single fleet.  Initial population size (N0) and catchability coefficients 
with standard errors are estimated from cpue and catch, respectively, via multiple linear regression (Leslie-Davis) 
and nonlinear regression (catch equation).
}
\usage{
depletM(catch = NULL, effort = NULL, M = NULL, 
 method = c("l", "c"),stq = NULL, stN0 = NULL)
}
\arguments{
  \item{catch}{the vector containing catches for each removal period (in sequential order).}
  \item{effort}{the vector containing effort associated with catch for each removal period. Rows must match those of catch.}
  \item{M}{the constant natural mortality rate}
  \item{method}{the depletion method: l=Leslie-Davis, c=catch equation. Default is both.}
  \item{stq}{initial guess for catchability if the catch equation method is used.}
  \item{stN0}{initial guess for N0 if the catch equation method is used.}
}
\details{
  For the Leslie-Davis method, cpue is calculated internally from catch and effort. The natural mortality rate must be scaled
to time period between removal (see Rosenberg et al. (1990) for more details). The \code{lm} and \code{nls} functions are the estimation  
methods for the Leslie-Davis and catch equations, respectively. 
}
\value{
  \item{results }{list element containing a table of parameters, estimates, associated standard errors, and 
 confidence intervals for each method}
  \item{outpt}{list element containing the calculated regressors, predicted cpue (Leslie) or catch (catch equation) and residuals}
}
\references{Rosenberg, A.A., G. P. Kirkwood, J. A. Crombie, and J. R. Beddington. 1990. The assessment of stocks of annual
squid species. Fish. Res. 335-350.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\keyword{misc}
