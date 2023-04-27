\name{growth_sel}
\alias{growth_sel}

\title{
Fitting a von Bertalanffy curve to length and age data biased by gear selectivity 
}
\description{
A von Bertalanffy growth curve is fitted to age and length data corrected for gear selectivity via the method of
Schueller et al. (2014). 
}
\usage{
growth_sel(age = NULL, size = NULL, weights = NULL, minlimit = NULL, maxlimit = NULL,
 minmax = NULL, switch_varpar = 1, 
 Linf = list(init = 1000, lb = 100, ub = 2000, prior.mean = 1000, prior.var = -0.5,
 prior.pdf = 1), 
 K = list(init = 0.3, lb = 0.1, ub = 0.9, prior.mean = 0.3, prior.var = -0.05,
 prior.pdf = 1), 
 t0 = list(init = -0.5, lb = -2, ub = -1e-04, prior.mean = -0.5, prior.var = -0.5,
 prior.pdf = 1), 
 varpar = list(init = 50, lb = 10, ub = 100, prior.mean = 5, prior.var = -1,
 prior.pdf = 1),
 tmb.control = list(maxit = 5000, trace = F),
 nlminb.control = list(eval.max = 1e+05, iter.max = 1000),
 species_info = list(species = NULL, size_units = NULL))
}
\arguments{
  \item{age}{a vector of ages.}
  \item{size}{a vector of body sizes associated with the age data.}
  \item{weights}{a vector of observation weights associated with length data and used to produce weighted likelihood. 
		Set to 1 for unweighted likelihood.}
  \item{minlimit}{a single value or vector associated with the length data. If a single value, a vector the length of the age vector is produced.}
  \item{maxlimit}{a single value or vector associated with the length data. If a single value, a vector the length of the age vector is produced.}
  \item{minmax}{a vector of 1 and 2s indicating whether the data row is being applied to the minimum (1) or maximum part (2) of the likelihood. 
   In general, the break between a 1 and 2 would be the age that has the fullest distribution of length (a well sampled age class where 
    no bias correction is expected).}
  \item{switch_varpar}{estimated variance parameter: 1 = standard deviation (sigma),  
   2 = CV (sigma / mean), 3 = variance to mean ratio (sigma^2/mean)}
  \item{Linf}{list specifying the initial starting value (\emph{init}) of L-infinity, the parameter's lower (\emph{lb}) and upper bounds (\emph{ul})
    for box constraints, prior mean (\emph{prior.mean}), prior variance (\emph{prior.variance}) and prior distribution (\emph{pdf}). 
     \emph{pdf}: 1 = prior not used, 2 = lognormal, 3 = normal, 4 = beta.}
  \item{K}{list specifying same arguments for \emph{K} as \emph{Linf}.}
  \item{t0}{list specifying same arguments for \emph{t0} as \emph{Linf}.}
  \item{varpar}{list specifying same arguments for the estimated variance parameter (\emph{varpar}) as \emph{Linf}.}
  \item{tmb.control}{controls for the \emph{MakeADFun} function. See package \emph{TMB} for more information.}
  \item{nlminb.control}{controls for the \emph{nlminb} function. See function \emph{nlminb} for more information.}
  \item{species_info}{list specifying the species analyzed (\emph{species}) and units of the size measurements (\emph{size_units}).}
}
\details{
 The von Bertalanffy growth model \emph{Lage=Linf*(1-exp(-K*(age-t0))} is fitted to length-at-age data 
  adjusted for bias related to selectivity of gears used to collect the length and age samples following the method of Schueller et al. (2014).
}
\value{
List containing list elements of the run information (\emph{run_info}), filtering indicator (\emph{message}), convergence information (\emph{convergence_info}),
 parameter estimates with associated standard errors and boundary values (\emph{estimates}), likelihood values (\emph{likelihood}) and 
 predicted values (\emph{predicted}).
}
\references{
 Schueller, A. M., E. H. Williams and R. T. Cheshire. 2014. A proposed, tested, and applied adjustment to account
for bias in growth parameter estimates due to selectivity. Fisheries Research 158: 26-39.
}
\author{
 Amy M. Schueller, National Marine Fisheries Service, Beaufort, NC \email{amy.schueller@noaa.gov}
}
\note{
Amy Schueller provided her AD Model Builder code which was translated to TMB code by Gary Nelson.
}
\examples{
 \dontrun{
  data(simulus)
  growth_sel(age=simulus$age,size=simulus$size,weights=simulus$weight,
    minlimit=simulus$minlimit,
    maxlimit=simulus$maxlimit,minmax=simulus$minmax,
    switch_varpar=1,
    Linf=list(init=1000,lb=100,ub=2000,prior.mean=1000,prior.var=-0.5,prior.pdf=1),
    K=list(init=0.3,lb=0.1,ub=0.9,prior.mean=0.3,prior.var=-0.05,prior.pdf=1),
    t0=list(init=-0.5,lb=-4,ub=-0.001,prior.mean=-0.5,prior.var=-0.5,prior.pdf=1),
    varpar=list(init=50.0,lb=10,ub=100,prior.mean=100,prior.var=-1.0,prior.pdf=1),
    tmb.control=list(maxit=5000,trace=F),nlminb.control=list(eval.max=100000,
    iter.max=1000),
    species_info=list(species="gag",size_units="inches"))
 }
}
\keyword{misc}
