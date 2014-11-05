\name{irm_h}
\alias{irm_h}
\title{Age-Independent Instantaneous Rates Tag Return Model of Hoenig et al. (1998)}
\description{
The age-independent instantaneous rates model of Hoenig et al. (1998) for estimating fishing and natural 
mortality from tag returns of harvested fish is implemented assuming known values of initial tag survival (phi) and 
reporting rate (lambda)
}
\usage{
irm_h(relyrs = NULL, recapyrs = NULL, N = NULL, recapharv = NULL,
lambda = NULL,phi = NULL, Fyr = NULL, Myr = NULL, initial = NULL,
lower = c(0.0001,0.0001),upper = c(5,5), maxiter = 10000)
}
\arguments{
  \item{relyrs}{vector containing the start and end year of the entire release period 
  (e.g., c(1992, 2006)).}
  \item{recapyrs}{vector containing the start year and end year of entire recapture period (e.g., c(1992, 2008)).}
  \item{N}{vector of total number of tagged fish released in each release year (one value per year).}
  \item{recapharv}{matrix of the number of tag recoveries of harvested fish by release year (row) and recovery year
  (column). The lower triangle (blank cells) may be filled with -1s as place holders.  Missing values in the upper 
   triangle (release/recovery cells) are not allowed.}
  \item{lambda}{vector of reporting rate estimates for harvested fish. One value for each recovery year.}
  \item{phi}{vector of initial tag survival estimates (phi) for harvested fish. One value for each recovery year.
    1=no loss}
  \item{Fyr}{vector of year values representing the beginning year of a period over which to estimate a constant 
   fishing mortality rate
 (F). If estimation of F for each recovery year is desired, enter the year value for each year. The first year value 
  must be the start year for the recovery period.}
  \item{Myr}{vector of year values representing the beginning year of a period over which to estimate a constant natural 
mortality rate (M). If estimation of M for each recovery year is desired, enter the year value for each year. The first 
year value must be the 	start year for the recovery period.}
  \item{initial}{vector of starting values for fishing, and	natural mortality estimates. 
First position is the starting value for all Fs and second position is the starting value for all Ms (e.g., c(0.1,0.2)).}
  \item{lower}{vector of lower bounds of F and M estimates used in optimization routine. 
First position is the lower value for all Fs and second position is the lower value for all Ms. Default = 0.0001.}
  \item{upper}{vector of upper bounds of F and M estimates used in optimization routine. 
First position is the upper value for all Fs and second position is the upper value for all Ms. Default = 5}
 \item{maxiter}{maximum number iterations used in the optimization routine.}
}
\details{
The instantaneous tag return model of Hoening et al. (1998) assuming known initial tag loss and reporting rates is
implemented. This model assumes that tagged fish are fully-recruited to the fishery and that fishing took place 
throughout the year.  The observed recovery matrices are compared to expected recovery matrices to estimate model 
parameters. Asymmetric recovery matrices are allowed (recovery years > release years). All summary statistics follow 
Burnham and Anderson (2002).  Model degrees of freedom are calculated as the number of non-zero cells from the 
harvested and released recovery matrices minus the number of estimated parameters. Total chi-square is calculated by 
summing cell chi-square values for all non-zero observed cells of the harvest, released, and not seen matrices.
C-hat, a measure of overdispersion, is estimated by dividing the total chi-square value by the model degrees of freedom.
Pooling of cells to achieve an expected cell value of 1 is performed and pooled chi-square and c-hat metrics are 
additionally calculated. Pearson residuals are calculated by subtracting the observed numbers of recoveries in each 
cell from the predicted numbers of recoveries and dividing each cell by the square-root of the predicted cell value. 
The \code{optim} routine is used to find the parameters that minimize the -1*negative
log-likelihood. The variance of instantaneous total mortality (Z) is calculated by \code{varF + varM + 2cov(F,M)}, and
the variance of survival (S) is estimated from the variance of Z using the delta method.

The program allows the configuration of different model structures (biological realistic models) for the estimation of
fishing and natural mortalities.  Consider the following examples:

\emph{Example 1}

Release years range from 1991 to 2003 and recovery years from 1991 to 2003.  One model structure might be constant 
fishing mortality estimates over the recovery years of 1991-1994 and 1995-2003, and one constant estimate of natural
mortality for the entire recovery period.  To specify this model structure, the beginning year of each interval is 
assigned to the \code{Fyr} vector (e.g.,Fyr<-c(1991, 1995)), and the beginning year of the recovery period 
is assigned to the \code{Myr} vector  (e.g.,Myr<-c(1991)).  The first value of each vector must always be
 the beginning year of the recovery period regardless of the model structure.

\emph{Example 2}

Release years range from 1991 to 2003 and recovery years from 1991 to 2003.  One model might be fishing mortality 
estimates for each year of recovery years and two constant estimates of natural mortality for 1991-1996 and
1997-2003.  To specify this model structure, one value for each year is assigned to the Fyr vector 
(e.g., Fyr<-c(1991,1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003) and the beginning 
years of the natural mortality intervals are assigned to the Myr vector (e.g.,Myr<-c(1991, 1997)).

Averaging of model results can be accomplished using the function \code{tag_model_avg}.
}
\value{List containing summary statistics for the model fit, model convergence status, parameter correlation matrix,
estimates of fishing mortality, natural mortality, total instantaneous mortality (Z), and survival (S)
and their variances and standard errors by year, observed and predicted recoveries for harvested, released, and
 "not-seen" fish, cell chi-square and Pearson values for harvested, released, and "not seen" fish and a model 
configuration label (type) used in the \code{tag_model_avg} function.

}

\references{ 
Burnham, K. P. and D. R. Anderson. 2002. Model selection and multimodel inference : A Practical Information-Theorectic
 Approach, 2nd edition. Spriner-Verlag, New York, NY. 488 p.

Hoenig, J. M, N. J. Barrowman, W. S. Hearn, and K. H. Pollock. 1998. Multiyear tagging studies incorporating fishing 
effort data.  Canadian Journal of Fisheries and Aquatic Sciences 55: 1466-1476.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{irm_cr}} \code{\link{tag_model_avg}}}
\examples{
# Data come from Table 4 and model structure from Table 5 under "year-specific F, 
# constant M" in Hoenig et al. (1998)  
data(Hoenig)
model1<-irm_h(relyrs = Hoenig$relyrs, recapyrs = Hoenig$recapyrs, 
N = Hoenig$N, recapharv = Hoenig$recapharv,lambda = Hoenig$lambda,
phi = Hoenig$phi, Fyr = Hoenig$Fyr, Myr = Hoenig$Myr, initial = c(0.1,0.1), 
lower = c(0.0001,0.0001),upper = c(5,5), maxiter = 10000)
}            
\keyword{misc}
