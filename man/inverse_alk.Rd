\name{inverse_alk}
\alias{inverse_alk}
\title{Inverse Age-Length Key Method of Hoenig and Heisey (1987)}
\description{
  Estimate the age composition of fish from size frequencies with missing age data by using the inverse method of Hoenig and Heisey (1986).
}
\usage{
inverse_alk(alk1 = NULL, lf1 = NULL, lf2 = NULL, toler = 0.000001, max.iter = 10000)
}
\arguments{
  \item{alk1}{an age-size matrix (numbers) with age as rows and size intervals as columns.}
  \item{lf1}{an optional vector of number of fish at size that will be used to expand the numbers in \emph{alk1}
     matrix to the numbers in the vector before estimation of age composition of lf2. 
     Vector length must match number of columns in \emph{alk1}. \emph{NULL} indicates no vector used.}
  \item{lf2}{a required vector of number of fish at size for which age composition will be estimated. Vector length must 
    match number of columns in \emph{alk1}.}
\item{toler}{convergence criterion. The iterations end when \emph{|L'-L|} <= \emph{toler}. Default = 0.000001.}
\item{max.iter}{additional convergence criterion. The maximum number of iterations allowed. Default =10000.}
}

\details{
 The inverse age-length key method of Hoenig and Heisey (1987) is used to estimate age composition of a sample of size data with no age data from
an age-length key which may be from a different year/ region. The method estimates the probability of size given age which
is not affected by variability in recruitment and survival (Ailloud and Hoenig, 2019). What does affect the probability of size given
age is spatiotemporal variations of size at age. These could be caused by changes in growth rates, or changes in mean size at age
due to changes in fishing practices, for example. So the inverse key can be applied to samples from populations with differing
age compositions than the population from which it was derived, so long as size at age does not vary considerably among sampling
events (copied from Ailloud and Hoenig, 2019).
  }
\value{list containing observed objects (alk1, lf1, lf2), the estimated alk in numbers for lf2, residuals, 
and the estimated age composition for lf2.}
\references{
Ailloud, L. E. and J. M. Hoenig. 2019. A general theory of age-length keys: combining the forward and inverse keys to estimate age composition from
incomplete data. ICES Journal of Marine Science. 76: 1515-1523.

Hoenig, J. M. and D. M. Heisey. 1987. Use of a log-linear model with the EM 
algorithm to correct estimates of stock composition and to convert length to age.
Transactions of the American Fisheries Society 116: 232-243.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{alkD}} \code{\link{alkss}} \code{\link{alkprop}}}
\examples{
 \dontrun{
  data(hohe) 
  inverse_alk(alk1=hohe$alk,lf1=NULL,lf2=hohe$lengths,toler=0.000001,max.iter=100000)
 }
}
\keyword{misc}
