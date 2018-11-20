\name{sr}
\alias{sr}
\title{Estimation and Model Comparison of Stock-Recruitment Relationships}
\description{
This function fits 14 models of recruitment-stock relationships
to recruitment numbers and spawning stock (e.g., spawning stock biomass or fecundity) data and provides 
model selection statistics for determining the best model fit.}
\usage{
sr(recruits = NULL, stock = NULL, model = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
 10, 11, 12, 13, 14), 
select = 1, initial = list(RA = NULL, RB = NULL, Rrho = NULL, BHA = NULL, 
BHB = NULL, BHrho = NULL,
 SHA = NULL, SHB = NULL, SHC = NULL, DSA = NULL, DSB = NULL, DSC = NULL, 
MYA = NULL, MYB = NULL, 
 MYC = NULL), control = list(maxit = 10000), plot = FALSE)
}
\arguments{
  \item{recruits}{a vector of numbers of recruits}
  \item{stock}{any spawning stock quantity (e.g., spawning biomass, numbers, fecundity) corresponding to the
   vector of recruits.}
  \item{model}{the model to fit. Models are
   0 = Density-Independent, 1 = Ricker with uncorrelated normal errors (N-U), 2 = Ricker with uncorrelated log-normal errors (L-U),
   3 = Ricker with correlated normal errors (N-C), 4 = Ricker with correlated log-normal errors (L-C),
   5 = Beverton-Holt with uncorrelated normal errors, 6 = Beverton-Holt with uncorrelated log-normal errors,
   7 = Beverton-Holt with correlated normal errors, 8 = Beverton-Holt with correlated log-normal errors,
   9 = Shepherd with uncorrelated normal errors, 10 = Shepherd with uncorrelated log-normal errors, 11 = Deriso-Schnute with
       uncorrelated normal errors, 12 = Deriso-Schnute with uncorrelated log-normal errors, 12 = Myers depensatory model with 
      uncorrelated normal errors, and 14 = Myers depensatory model with uncorrelated log-normal errors. Default is all.
  }
  \item{select}{method used to determine starting values. 1 = automatic, 2 = user-specified. Default=1. Automatic
	selection of starting might not always work given the data provided.}
  \item{initial}{if select = 2, list of starting values for each equation type. See equation parameter names in \emph{Details}.}
  \item{control}{see function \emph{optim}.}
  \item{plot}{logical indicating whether an observed-predicted plot should be produced. Default = FALSE.}
}
\details{

The following equations are fitted:

Ricker: \code{recruits = RA*stock*exp(-RB*stock)}

Beverton-Holt: \code{recruits = (BHA*stock)/(1+(BHA*stock)/BHB)}

Shepherd: \code{recruits = (SHA*stock)/(1+SHB*stock^SHC)}

Deriso-Schnute: \code{recruits = DSA*stock*(1-DSB*DSC*stock)^(1/DSC)}

Myers: \code{(MYA*datar$stock^MYC)/(1+((datar$stock^MYC)/MYB))} 

Maximum likelihood is used to estimate model parameters. 

For uncorrelated normal errors, the negative log-likelihood is

\code{n/2*log(2*pi)+n*log(sqrt(sigma2))+1/(2*sigma2)*sum((recruits-predicted)^2)}

where n is the number of observation, sigma2 is the maximum likelihood of residual variance and
predicted is the model predicted recruits.  sigma2 is calculated internally as 

\code{sigma2 = sum((recruits-predicted)^2)/n}.


For uncorrelated log-normal errors, the negative log-likeliood is

\code{n/2*log(2*pi)+n*log(sqrt(lsigma2))+sum(log(recruits))+1/(2*lsigma2)*}


\code{sum((log(recruits)-log(predicted)+lsigma2/2)^2)}

lsigma2 is calculated internally as \code{lsigma2 = sum((log(recruits)-log(predicted))^2)/n}.

For correlated normal errors, the negative log-likelihood is

\code{n/2*log(2*pi)+n*log(sqrt(sigma2w))-0.5*log(1-rho^2)+} 

\code{1/(2*sigma2w)*sumR+((1-rho^2)/(2*sigma2w))*(datar$recruits[1]-predicted[1])^2}

where rho is the estimated autocorrelation (AR1) parameter, sigma2w is the white noise residual variance, and sumR is calculated as

\code{for(k in 2:n) sumR<-sumR+(recruits[k]-rho*recruits[k-1]-}

\code{predicted[k]+rho*predicted[k-1])^2}

sigma2w is calculated internally as

    \code{res = recruits - predicted}

     \code{es = c(res[1:c(length(res)-1)]*rho)}

     \code{sigma2w = sum((res[-1]-es)^2)/c(n-1)}
 

For correlated log-normal errors, the negative log-likelihood is 

 \code{n/2*log(2*pi)+n*log(sqrt(lsigma2w))+sum(log(recruits))-0.5*log(1-rho^2)+}

\code{1/(2*lsigma2w)*lsumR+((1-rho^2)/(2*lsigma2w))*(log(recruits[1])-}

\code{log(predicted[1])+lsigma2w/2)^2}
     
where lsumR is calculated as
     
\code{for(k in 2:n) lsumR<-lsumR+(log(recruits[k])-pho*log(recruits[k-1])}

\code{-log(predicted[k])+rho*log(predicted[k-1])+(1-phi)*lsigma2w/2)^2}

and lsigma2w is calculated as

  \code{res = log(recruits)-log(predicted)}

  \code{es = c(res[1:c(length(res)-1)]*pho)}

  \code{lsigma2w = sum((res[-1]-es)^2)/c(n-1)}.

Correlated error structures are available for the Ricker and Beverton-Holt model only.  The
 names for specification of starting values of the AR1 parameter are \code{Rrho} and \code{BHrho}.  

Akaike Information Criterion for small sample sizes (AICc), Akaike weights and evidence ratios (Burham and Anderson 2002)
 are provided for each model selected above. 

This function uses function \emph{optim} to estimate parameters and function \emph{hessian} in package \emph{numDeriv} to calculate the
hessian matrix from which standard errors are derived.
}
\value{
 Lists containing estimation results. \emph{results} contains parameter estimates,
 associated standard errors, residual variances, negative log-likelihoods and AICc values for each model.
 If the standard errors are \code{NaN}, the hessian could not be inverted (i.e., poor model fit). 
\emph{evidence_ratios} contains Akaike weights and evidence ratios for model selection.
\emph{convergence} contains convergence criterion: 0 = no problems, >0 = problems (see function \emph{optim}).
\emph{correlations} contains the estimated parameter correlations. Correlation will be NA if hessian could not be
inverted. \emph{predicted} contains the predicted values from each model. \emph{residuals} contains the residuals from each model.
}
\references{

Brodziak, J, and C. M. Legault. 2005. Model averaging to estimate rebuilding targets for
overfished stocks. Canadian Journal of Fisheries and Aquatic Sciences 62: 544-562.

Brodziak, J, and C. M. Legault. 2010. Reference manual for SRFIT version 7. NOAA Fisheries Toolbox.

Burnham, K. P. and D. R. Anderson. 2002. Model Selection and Multimodel Inference, Second edition.
Springer-Verlag New York, New York.  488 pages.
 
Myers, R. A., N. J. Barrowman, J. A. Hutching and A. A. Rosenberg. 1995. Population dynamics
of exploited fish stocks at low population levels. Science 269: 1106-1108.

Quinn, T. J. and R. B. Deriso. 1999. Quantitative fish dynamics. Oxford University Press. 542 pages. 

}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}

\examples{
 \dontrun{
data(striper)
outs<-sr(recruits=striper$recruits,stock=striper$stock,select=2,model=c(5,6,7,8),
         initial=list(RA=5e3,RB=2e-5,Rrho=0.1,
                      BHA=8e3,BHB=1e8,BHrho=0.1,
                      SHA=1.5e3,SHB=5.6e8,SHC=1,
                      DSA=9e3,DSB=9e-5,DSC=-1.14,
                      MYA=1e6,MYB=1e5,MYC=0.4),plot=TRUE)
}
}
\keyword{misc}
