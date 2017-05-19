\name{mature}
\alias{mature}
\title{Estimation of proportion mature at age when immature fish are unavailable}
\description{
  Calculates proportion mature-at-age based on Maki et al. (2001).
}
\usage{
mature(cap_age=NULL, mature_age=NULL, age_all_immature=NULL,
age_all_mature=NULL, initial=NULL, nrandoms=1000)
}
\arguments{
  \item{cap_age}{vector of ages representing age when fish was capture. 
 	One record per individual.}
  \item{mature_age}{vector of ages representing age at which individual 
   mature.One record per individual.}
  \item{age_all_immature}{age at which all fish are deemed immature. All ages below this age 
  are assumed immature also.}
  \item{age_all_mature}{age at which all fish are deemed mature.
	All ages above this age are also assumed mature.}
  \item{initial}{starting values for proportion estimates. There should be
\emph{age_all_mature - age_all_immature-2} values. If not, the last value is used for 
missing values or if the vector is too large, the vector is truncated.}
  \item{nrandoms}{the number of randomizations used to estimate standard errors.}
}
\details{
 Estimation of probability follows Maki et al. (2001).The standard errors of
 parameters are estimated via Monte Carlos methods where the number of each maturing age 
 for each capture age are randomly draw from a multinomial distribution parameterized with 
 probabilities and total sample size of the original data. The methods of Maki et al. (2001)
 are applied to the randomized data and the randomization is repeated \emph{nrandoms} times.
  The mean and standard deviation of all runs are treated as the parameter estimates and standard errors. 
}
\value{a list object containing the estimated proportions-at-age and standard errors,
     the original data and expected values
}
\references{Maki, K. L., J. M. Hoenig and J. E. Olney. 2001. Estimating proportion mature at age 
when immature fish are unavailable for study, with applications to American shad in the York
River, Virginia. North Am. J. Fish. Manage. 21: 703-716.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
  \dontrun{
   ## Maki data for 3 years combined
   data(maki)
    mature(cap_age=maki$capture_age,mature_age=maki$age_mature,age_all_immature=2,
                 age_all_mature=8,initial=c(0.1,0.05,0.05,0.05),nrandoms=1000)
  }

}
\keyword{misc}

