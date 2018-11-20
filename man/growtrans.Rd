\name{growtrans}
\alias{growtrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Growth Transition Matrix for a Size-Structured Population Dynamics Model }
\description{
Generates a growth transition matrix from parameters of the von Bertalanffy growth 
equation following Chen et al. (2003)}
\usage{
growtrans(Lmin = NULL, Lmax = NULL, Linc = NULL, Linf = NULL, SELinf = NULL,
          K = NULL, SEK = NULL, rhoLinfK = NULL)
}
\arguments{
  \item{Lmin}{Mid-point of starting size class.}
  \item{Lmax}{Mid-point of end size class.  This should be one increment larger than Linf.}
  \item{Linc}{Size class increment.}
  \item{Linf}{L-infinity parameter of the von Bertalanffy growth equation.}
  \item{SELinf}{Standard error of Linf.}
  \item{K}{Growth parameter of the von Bertalanffy growth equation.}
  \item{SEK}{Standard error of K.}
  \item{rhoLinfK}{Correlation between Linf and K. Usually from a parameter correlation matrix.}
}
\details{
  Transition probabilities are calculated by using formulae 3-9 and procedures 
  in Chen et al. (2003). Negative growth increments result if Lmax is beyond 
  Linf, so the transition matrix is truncated at Linf. The last size class acts as 
  a plus group and has a probability of 1.
}
\value{A matrix of dimensions n size classes x n size classes.}
\references{
Chen, Y., M. Hunter, R. Vadas, and B. Beal. 2003. Developing a growth-transition matrix for stock assessment 
of the green sea urchin (Strongylocentrotus droebachiensis) off Maine. Fish. Bull. 101: 737-744.
}
\note{This function is based on an example EXCEL spreadsheet provided by Yong Chen.}   
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\examples{
# For Chen et al. 2003
growtrans(Lmin=40,Lmax=101,Linc=1,Linf=100,SELinf=15,K=0.100588,SEK=0.04255,rhoLinfK=0.94)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
