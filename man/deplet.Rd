\name{deplet}
\alias{deplet}
\title{Catch-Effort Depletion Methods For a Closed Population}
\description{
  Variable and constant effort models for the estimation of abundance from catch-effort depletion data assuming a closed population.
}
\usage{
deplet(catch = NULL, effort = NULL, method = c("l", "d", "ml",
 "hosc", "hesc", "hemqle", "sch"), nboot = 500)
}
\arguments{
\item{catch}{the vector containing catches for each removal period (in sequential order).}
\item{effort}{the vector containing effort associated with catch for each removal period. Rows must match those of catch.}
\item{method}{the depletion method. \emph{Variable Effort Models}: \code{l}= Leslie estimator, \code{d}= effort corrected
 Delury estimator, \code{ml}= maximum likelihood estimator of Gould and Pollock (1997), \code{hosc}= sampling coverage
 estimator for homogeneous model of Chao and Chang (1999), \code{hesc}= sampling coverage estimator for heterogeneous 
model of Chao and Chang (1999), and \code{hemqle}= maximum quasi likelihood estimator for heterogeneous model of
 Chao and Chang (1999). \emph{Constant Effort Models}: \code{sch}= maximum likelihood models that test 
for constant catchability.
}
\item{nboot}{the number of bootstrap resamples for estimation of standard errors in the \code{ml},
\code{hosc},\code{hesc}, and \code{hemqle} methods}
}
\details{
  The variable effort models include the Leslie-Davis (\code{l}) estimator (Leslie and Davis, 1939), the effort-corrected Delury (\code{d}) estimator (Delury,1947; Braaten, 1969),
the maximum likelihood (\code{ml}) method of Gould and Pollock (1997), sample coverage estimator for the homogeneous model (\code{hosc}) of Chao and Chang (1999), 
sample coverage estimator for the heterogeneous model (\code{hesc}) of Chao and Chang (1999), and the maximum quasi-likelihood estimator for the heterogeneous model (\code{hemqle}) of Chao and Chang (1999). 
The variable effort models can be applied to constant effort data by simply filling the \code{effort} vector with 1s.
The constant effort models include only model 1 (constant catchability) and model 2 (different catchability in first period) (\code{sch}) of Schnute (1983) which are equivalent to the generalized removal models for k=1 and k=2, respectively,
of White, Anderson, Burnham, and Otis (1982:p. 111-114).  A vector of effort data is not required for the Schnute models.  Note: Calculation of the standard error using the \code{ml} method takes considerable time.
}
\value{
  Separate output lists with the method name and extension \code{.out} are created for each method and contain tables of various statistics associated with the method.  
}
\references{Braaten, D. O. 1969. Robustness of the Delury population estimator. J. Fish. Res. Board Can. 26: 339-355.

Chao, A. and S. Chang. 1999. An estimating function approach to the inference of catch-effort models. Environ. Ecol. Stat. 6: 313-334.

Delury, D. B. 1947. On the estimation of biological populations. Biometrics 3: 145-167.

Gould, W. R. and H. H. Pollock. 1997. Catch-effort maximum likelihood estimation of important population parameters. Can. J. Fish. Aquat. Sci 54: 890-897.

Leslie, P. H. and D. H.S. Davis. 1939. An attempt to determine the absolute number of rats on a given area.  J. Anim. Ecol. 9: 94-113.

Schnute, J. 1983. A new approach to estimating populations by the removal method. Can. J. Fish. Aquat. Sci. 40: 2153-2169.

White, G. C., D. R. Anderson, K. P. Burnham, and D. L. Otis. 1982. Capture-recapture and Removal Methods for Sampling
Closed Populations. Los Alamos National Laboratory LA-8787-NERP. 235 p.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
data(darter)
deplet(catch=darter$catch,effort=darter$effort,method="hosc") 
hosc.out
}            
\keyword{misc}
