\name{catchsurvey}
\alias{catchsurvey}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Catch-Survey Analysis}
\description{
 This function applies the catch-survey analysis method of Collie and Kruse (1998) for estimating abundance
 from catch and survey indices of relative abundance 
}
\usage{
catchsurvey(year = NULL, catch = NULL, recr = NULL, post = NULL, M = NULL,
 T = NULL, phi = NULL, w = 1, initial = c(NA,NA,NA),uprn = NA, graph = TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{vector containing the time series of numeric year labels.}
  \item{catch}{vector containing the time series of catch (landings) data.}
  \item{recr}{vector containing the time series of survey indices for recruit individuals.}
  \item{post}{vector containing the time series of survey indices for post-recruit individuals.}
  \item{M}{instantaneous natural mortality rate. Assumed constant throughout time series}
  \item{T}{proportion of year between survey and fishery.}
  \item{phi}{relative recruit catchability.}
  \item{w}{recruit sum of squares weight.}
  \item{initial}{initial recruit estimate,initial postrecruit estimate in year 1, and
         initial catchability estimate.}
  \item{uprn}{the upper bound for the recruit and postrecruit estimates.}
  \item{graph}{logical indicating whether observed versus predicted recruit and post-recruit indices,
    total abundance and fishing mortality should be plotted. Default=TRUE.}
}
\details{
 Details of the model are given in Collie and Kruse (1998).    
}
\value{List containing the estimate of catchability, predicted recruit index by year (rest), estimate of recruit abundance (R),
 predicted post-recruit index by year (nest), post-recruit abundance (N), total abundance (TA: R+N), total instantaneous mortality (Z),
 and fishing mortality (Fmort)}
\references{Collie JS and GH Kruse 1998. Estimating king crab (Paralithodes camtschaticus) abundance from
commercial catch and research survey data. In: Jamieson GS, Campbell A, eds. Proceedings of the North Pacific Symposium 
on Invertebrate Stock Assessment and Management. Can Spec Publ Fish Aquat Sci. 125; p 73-83.

See also Collie JS and MP Sissenwine. 1983. Estimating population size from relative abundance data
measured with error. Can J Fish Aquat Sci. 40:1871-1879.
}
\author{ Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
 ## Example take a bit of time to run
  \dontrun{
   data(nshrimp)
   catchsurvey(year=nshrimp$year,catch=nshrimp$C,recr=nshrimp$r,post=nshrimp$n,M=0.25,T=0.5,phi=0.9,w=1,
       initial=c(500,500,0.7),uprn=10000)}
   }
\keyword{misc}

