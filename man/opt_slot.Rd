\name{opt_slot}
\alias{opt_slot}
\title{Optimum Slot and Trophy Size Limits for Recreational Fisheries}
\description{
    Calculates optimum trophy catch given a slot size over a range of F values. Also, finds Fmax for a cohort given age-at-first recruitment,
age-at-first-entry, slot age, and age at which fish are considered trophy size following Jensen (1981).

}
\usage{
opt_slot(M = NULL, N = 1000, recage = NULL, entage = NULL,
 trage = NULL, slage = NULL,  stF = 0, endF = 2, intF = 0.05)
}
\arguments{
 \item{M}{natural mortality}
  \item{N}{ cohort size}
  \item{recage}{age-at-first recruitment}
  \item{entage}{age-at-entry into the fishery}
  \item{slage}{upper age of slot for legal fish}
  \item{trage}{age of fish considered trophy size}
  \item{stF}{starting F of range to explore}
  \item{endF}{ending F of range to explore}
  \item{intF}{increment of F}
}
\details{
   Calculations follow equations given in Jensen (1981).
}
\value{
 \item{Catch}{dataframe containing range of Fs and associated total catch, nontrophy, and trophy catch of designated cohort size}
  \item{Fmax }{F at which trophy catch is maximum given slot}
}
\references{Jense, A. L. 1981. Optimum size limits for trout fisheries. Can. J. Fish. Aquat. Sci. 38: 657-661.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{opt_trophy}}}
\examples{
# Example from Jensen (1981) page 661
opt_slot(M=0.70,N=1000,recage=1,entage=1,slage=3,trage=4)
}
\keyword{misc}