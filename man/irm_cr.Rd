\name{irm_cr}
\alias{irm_cr}
\title{Age-Independent Instantaneous Rates Model of Jiang et al. (2007) Incorporating Catch and Release Tag Returns }
\description{
The age-independent instantaneous rates model of Jiang et al. (2007) for estimating fishing and natural 
mortality from catch-release tag returns is implemented assuming known values of initial tag survival (phi) and 
reporting rate (lambda)
}
\usage{
irm_cr(relyrs = NULL, recapyrs = NULL, N = NULL, recapharv = NULL, 
recaprel = NULL, hlambda = NULL, rlambda = NULL, hphi = NULL, 
rphi = NULL, hmrate = NULL, Fyr = NULL, FAyr = NULL, Myr = NULL,
initial = c(0.1,0.05,0.1), lower = c(0.0001,0.0001,0.0001), 
upper=c(5,5,5),maxiter=10000)
}
\arguments{
  \item{relyrs}{vector containing the start and end year of the entire release period (e.g., c(1992,2006)).}
  \item{recapyrs}{vector containing the start year and end year of entire recapture period (e.g., c(1992,2008)).}
  \item{N}{vector of total number of tagged fish released in each release year (one value per year).}
  \item{recapharv}{matrix of the number of tag recoveries of harvested fish by release year (row) and recovery year 
     (column). The lower triangle (blank cells) may be filled with -1s as place holders.  Missing values in the upper
     triangle (release/recovery cells) are not allowed.}
  \item{recaprel}{matrix of the number of tag recoveries of fish recaptured and re-released with the tag removed 
      by release year (row) and recovery year (column). The lower triangle (blank cells) may be filled with -1s as 
      place holders.  Missing values in the upper triangle (release/recovery cells) are not allowed.}
  \item{hlambda}{vector of reporting rate estimates (lambda) for harvested fish. One value for each recovery year.}
  \item{rlambda}{vector of reporting rate estimates (lambda) for recaptured fish re-released with tag removed. 
    One value for each recovery year.}
  \item{hphi}{vector of initial tag survival estimates (phi) for harvested fish. One value for each recovery year.
   1 = no loss}
  \item{rphi}{vector of initial tag survival estimates (phi) for recaptured fish re-released with tag removed fish. 
     One value for each recovery year. 1 = no loss}
  \item{hmrate}{vector of hooking mortality rates. One value for each recovery year.}
  \item{Fyr}{vector of year values representing the beginning year of a period over which to estimate a constant 
     fishing mortality rate (F). If estimation of F for each recovery year is desired, enter the year value for each
     year. The first year value must be the start year for the recovery period.}
  \item{FAyr}{vector of year values representing the beginning year of a period over which to estimate a constant tag 
     mortality rate (FA). If estimation of FA for each recovery year is desired, enter the year value for each year. 
     The first year value must be the start year for the recovery period.}
  \item{Myr}{vector of year values representing the beginning year of a period over which to estimate a constant natural 
mortality rate (M). If estimation of M for each recovery year is desired, enter the year value for each year. The first 
year value must be the 	start year for the recovery period.}
  \item{initial}{vector of starting values for fishing, tag, and 	natural mortality estimates. First position is the 
   starting value for all Fs, second position is the starting value for all FAs, and the third position is the starting 
   value for all Ms (e.g., c(0.1,0.1,0.2)).}
  \item{lower}{vector of lower bounds of F, FA, and M estimates used in optimization routine. First position is the 
   lower value for all Fs, second position is the lower value for all FAs, and the third position
    is the lower value for all Ms.}
  \item{upper}{vector of upper bounds of F, FA, and M estimates used in optimization routine. First position is the 
   upper value for all Fs, second position is the upper value for all FAs, and the third position 
   is the upper value for all Ms.}
 \item{maxiter}{maximum number iterations used in the optimization routine.}
}
\details{
Jiang et al (2007) provides an extension of the Hoenig et al. (1998) instantaneous tag return model to account for 
catch/release of tagged fish. The benefits of this instantaneous rates model are that data from tagged fish that 
are recaptured and released alive are directly incorporated in the estimation of fishing and natural mortality.  
Jiang et al. models mortality of harvested fish and the mortality experienced by the tag because fish are often 
released after the tag has been removed. Therefore, additional tag mortality parameters are estimated in the model. 
The age-independent model of Jiang et al. is implemented here and initial tag loss and reporting rates are assumed 
known.  This model assumes that tagged fish are fully-recruited to the fishery and that fishing took place throughout
the year.  Similar to Hoenig et al. (1998), observed recovery matrices from the harvest and catch/release fish with 
removed tags are compared to expected recovery matrices to estimate model parameters. Asymmetric recovery matrices 
are allowed (recovery years > release years). All summary statistics follow Burnham and Anderson (2002).  Model 
degrees of freedom are calculated as the number of non-zero cells from the harvested and released recovery matrices 
minus the number of estimated parameters. Total chi-square is calculated by summing cell chi-square values for all 
non-zero observed cells of the harvest, released, and not seen matrices. C-hat, a measure of overdispersion, is 
estimated by dividing the total chi-square value by the model degrees of freedom. Pooling of cells to achieve an 
expected cell value of 1 is performed and pooled chi-square and c-hat metrics are additionally calculated.
Pearson residuals are calculated by subtracting the observed numbers of recoveries in each cell from the predicted 
numbers of recoveries and dividing each cell by the square-root of the predicted cell value. The variance of 
instantaneous total mortality (Z) is calculated by \code{varF + hmrate^2*varFA + varM + 2*sum(cov(F,M)+
hmrate^2*cov(F,FA)+hmrate^2*cov(FA,M))}, and the variance of survival (S) is calculated from Z using the delta method.  
The \code{optim} routine is used to find the parameters that minimize the -1*negative log-likelihood.

The program allows the configuration of different model structures (biological realistic models) for the estimation of
fishing, natural, and tag mortalities.  Consider the following examples:

\emph{Example 1}

Release years range from 1991 to 2003 and recovery years from 1991 to 2003.  One model structure might be constant 
fishing mortality estimates over the recovery years of 1991-1994 and 1995-2003, one constant estimate of tag mortality 
and one constant estimate of natural mortality for the entire recovery period.  To designate this model structure, 
the beginning year of each interval is assigned to the \code{Fyr} vector (e.g.,\code{Fyr<-c(1991, 1995)}), and the 
beginning year of the recovery period is assigned to the \code{FAyr} vector and the \code{Myr} vector 
(e.g., \code{FAyr<-c(1991)}; \code{Myr<-c(1991)}).  The first value of each vector must always be the beginning year
of the recovery period regardless of the model structure.

\emph{Example 2}

Release years range from 1991 to 2003 and recovery years from 1991 to 2003.  One model might be fishing and tag 
mortality estimates for each year of recovery years and two constant estimates of natural mortality for 1991-1996 and
1997-2003.  To designate this model structure, one value for each year is assigned to the Fyr and FAyr 
vectors (e.g., Fyr<-c(1991,1992,1993,1994,1995,1996,1997, 1998,1999,2000,2001,2002,2003 and
FAyr<-c(1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003)), and the 
beginning years of the natural mortality intervals are assigned to the Myr vector (e.g.,Myr<-c(1991,1997)).

Averaging of model results can be accomplished using the function \code{tag_model_avg}.
}
\value{List containing summary statistics for the model fit, model convergence status, parameter correlation matrix,
estimates of fishing mortality, natural mortality, tag mortality, total instantaneous mortality (Z), and survival (S)
and their variances and standard errors by year, observed and predicted recoveries for harvested, released, and
 "not-seen" fish, cell chi-square and Pearson values for harvested, released, and "not seen" fish, and a model 
configuration label (type) used in the \code{tag_model_avg} function.
}

\references{ 
Burnham, K. P. and D. R. Anderson. 2002. Model selection and multimodel inference : A Practical Information-Theorectic 
Approach, 2nd edition. Spriner-Verlag, New York, NY. 488 p.

Hoenig, J. M, N. J. Barrowman, W. S. Hearn, and K. H. Pollock. 1998. Multiyear tagging studies incorporating 
fishing effort data.  Canadian Journal of Fisheries and Aquatic Sciences 55: 1466-1476.

Jiang, H. 2005. Age-dependent tag return models for estimating fishing mortality, natural mortality and selectivity. 
Doctoral dissertation. North Carolina State  University, Raleigh.

Jiang, H., K. H. Pollock, C. Brownie, J. M. Hoenig, R. J. Latour, B. K. Wells, and J. E. Hightower. 2007.  Tag 
return models allowing for harvest and catch and release: evidence of environmental and management impacts on 
striped bass fishing and natural mortality rates.  North Amercian Journal of Fisheries Management 27:387-396.

}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{irm_h}} \code{\link{tag_model_avg}}}
\examples{
## Data come from Appendix Table A2 and model structure from model (a) in
## Table 3.2 of Jiang (2005) 
## Example takes a bit of time to run
  \dontrun{
  data(Jiang)
   model1<-irm_cr(relyrs = Jiang$relyrs, recapyrs = Jiang$recapyrs, 
     N = Jiang$N, recapharv = Jiang$recapharv, recaprel = Jiang$recaprel,
     hlambda = Jiang$hlambda, rlambda = Jiang$rlambda, hphi = Jiang$hphi,
     rphi = Jiang$rphi, hmrate = Jiang$hmrate, Fyr = Jiang$Fyr,
     FAyr = Jiang$FAyr, Myr = Jiang$Myr, initial = c(0.1,0.05,0.1), 
     lower = c(0.0001,0.0001,0.0001), upper=c(5,5,5),maxiter=10000)
  }
}            
\keyword{misc}
