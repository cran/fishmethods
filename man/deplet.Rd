\name{deplet}
\alias{deplet}
\title{Catch-Effort Depletion Methods For a Closed Population}
\description{
  Variable and constant effort models for the estimation of abundance from catch-effort depletion data assuming a closed population.
}
\usage{
deplet(catch = NULL, effort = NULL, method = c("l", "d", "ml",
 "hosc", "hesc", "hemqle", "wh"), kwh=NULL, nboot = 500, Nstart=NULL)
}
\arguments{
\item{catch}{the vector containing catches for each removal period (in sequential order).}
\item{effort}{the vector containing effort associated with catch for each removal period. Rows must match those of catch.}
\item{method}{the depletion method. \emph{Variable Effort Models}: \code{l}= Leslie estimator, \code{d}= effort corrected
 Delury estimator, \code{ml}= maximum likelihood estimator of Gould and Pollock (1997), \code{hosc}= sampling coverage
 estimator for homogeneous model of Chao and Chang (1999), \code{hesc}= sampling coverage estimator for heterogeneous 
model of Chao and Chang (1999), and \code{hemqle}= maximum quasi likelihood estimator for heterogeneous model of
 Chao and Chang (1999). \emph{Constant Effort Model}: \code{wh}= the generalized removal method of Otis et al. (1978).
}
\item{kwh}{the number of capture parameters (p) to fit in method \code{wh}. NULL for all possible capture parameters.}
\item{nboot}{the number of bootstrap resamples for estimation of standard errors in the \code{ml},
\code{hosc},\code{hesc}, and \code{hemqle} methods}
\item{Nstart}{starting value for N in method "wh". If NULL, start value is
automatically determined}
}
\details{
The variable effort models include the Leslie-Davis (\code{l}) estimator (Leslie and Davis, 1939), the effort-corrected Delury (\code{d}) estimator (Delury,1947; Braaten, 1969),
the maximum likelihood (\code{ml}) method of Gould and Pollock (1997), sample coverage estimator for the homogeneous model (\code{hosc}) of Chao and Chang (1999), 
sample coverage estimator for the heterogeneous model (\code{hesc}) of Chao and Chang (1999), and the maximum quasi-likelihood estimator for the heterogeneous model (\code{hemqle}) of Chao and Chang (1999). 
The variable effort models can be applied to constant effort data by simply filling the \code{effort} vector with 1s. Three removals are required to use the Leslie, Delury, and Gould
and Pollock methods.  

The constant effort model is the generalized removal method of Otis et al. 1978 reviewed in White et al. (1982: 109-114).
If only two removals, the two-pass estimator of N in White et al. (1982:105) and the variance estimator of Otis et al. (1978: 108) are used.

Note: Calculation of the standard error using the \code{ml} method may take considerable time.

For the Delury method, zero catch values are not allowed because the log-transform is used.

For the generalized removal models, if standard errors appear as \code{NA}s but parameter estimates are provided, the inversion of the Hessian failed.
If parameter estimates and standard errors appear as \code{NA}s, then model fitting failed.

For the Chao and Chang models, if the last catch value is zero, it is deleted from the data.  Zero values between positive values are permitted.
}
\value{
  Separate output lists with the method name and extension \code{.out} are created for each method and contain tables of various statistics associated with the method.  
}
\references{Braaten, D. O. 1969. Robustness of the Delury population estimator. J. Fish. Res. Board Can. 26: 339-355.

Chao, A. and S. Chang. 1999. An estimating function approach to the inference of catch-effort models. Environ. Ecol. Stat. 6: 313-334.

Delury, D. B. 1947. On the estimation of biological populations. Biometrics 3: 145-167.

Gould, W. R. and K. H. Pollock. 1997. Catch-effort maximum likelihood estimation of important population parameters. Can. J. Fish. Aquat. Sci 54: 890-897.

Leslie, P. H. and D. H.S. Davis. 1939. An attempt to determine the absolute number of rats on a given area.  J. Anim. Ecol. 9: 94-113.

Otis, D. L., K. P. Burnham, G. C. White, and D. R. Anderson. 1978. Statistical inference from capture data on closed animal populations. Wildl. Monogr. 62: 1-135.

White, G. C., D. R. Anderson, K. P. Burnham, and D. L. Otis. 1982. Capture-recapture and Removal Methods for Sampling
Closed Populations. Los Alamos National Laboratory LA-8787-NERP. 235 p.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
data(darter)
deplet(catch=darter$catch,effort=darter$effort,method="hosc") 
hosc.out
}            
\keyword{misc}
