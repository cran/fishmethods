\name{schnabel}
\alias{schnabel}
\title{Population Size Estimates from Repeated Mark-Recapture Experiments}
\description{
Estimates of population abundance from Schnabel (1938) and Schumacher and Eschmeyer (1943) are calculated from
repeated mark-recapture experiments following Krebs (1989). 
}
\usage{
schnabel(catch = NULL, recaps = NULL, newmarks = NULL,
 alpha = 0.05)
}
\arguments{
\item{catch}{A vector containing the number of animal caught in each mark-recapture experiment.}
\item{recaps}{A vector containing the number of animal recaptured in each mark-recapture experiment.}
\item{newmarks}{A vector containing the newly marked animals in each mark-recapture experiment.}
\item{alpha}{the alpha level for confidence intervals. Default = 0.05}
}
\details{
  All computations follow Krebs (1989: p. 30-34). For the Schnabel method, the poisson
distribution is used to set confidence intervals if the sum of all recaptures is <50,and the t distribution is used if the sum of all recaptures is >=50. 
For the Schumacher-Eschmeyer method, the t distribution is used to set confidence intervals.
}
\value{Dataframe containing the population estimates for the Schnabel and Schumacher & Eschmeyer methods (N),
the inverse standard errors (invSE), lower (LCI) and upper (UCI) confidence intervals,
 and the type of distribution used to set confidence intervals (CI Distribution). }
\references{ Krebs, C. J. 1989. \emph{Ecological Methodologies}. Harper and Row, New York, NY. 654 p.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
data(Gerking)
schnabel(catch=Gerking$C,recaps=Gerking$R, newmarks=Gerking$nM,
 alpha=0.10)
}
\keyword{misc}
