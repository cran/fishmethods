\name{compare2}
\alias{compare2}
\title{Comparisons of two age readers or two aging methods}
\description{
  Function compares graphically the readings of two age readers and calculates 2 chi-square statistics for 
  tests of symmetry.
}
\usage{
 compare2(readings, usecols = c(1,2), twovsone = TRUE, plot.summary = TRUE, 
  barplot = TRUE, chi = TRUE, pool.criterion = 1, cont.cor = TRUE, 
  correct = "Yates", first.name = "Reader A",second.name = "Reader B")
}
\arguments{
  \item{readings}{dataframe or matrix containing the readings by Reader 1 and those by Reader 2.}
  \item{usecols}{columns of the dataframe or matrix corresponding to the readings of Reader 1 and those of Reader 2.
	Default=c(1,2).}
  \item{twovsone}{logical for whether first type of graph is produced.}
  \item{plot.summary}{logical for whether summary table is put on first graph.}
  \item{barplot}{logical for whether barplot of frequency of disagreements is drawn.}
  \item{chi}{logical for whether 2 chi-square tests are performed.}
  \item{pool.criterion}{used to collapse pairs where the expected number of observations is < pooling
       criterion (default is 1).}
  \item{cont.cor }{logical for whether a continuity correction should be used in 1st chisquare test.}
  \item{correct}{character for whether "Yates" or "Edwards" continuity correction should be done (if cont.cor=TRUE).}
  \item{first.name}{character string describing the first reader or the first aging method. The
                 default is to specify "Reader A".}
  \item{second.name}{character string describing the second reader or the second aging method. The
                default is to specify "Reader B".}
}
\details{
   This function can plot the number of readings of age j by reader 2 versus the number of readings of age i by 
   reader 1 (if twovsone=TRUE). Optionally, it will add the number of readings above, on, and below the 45 degree 
   line (if plot.summary=TRUE). The function can make a histogram of the differences in readings (if barplot=TRUE). 
   Finally, the program can calculate 2 chi-square test statistics for tests of the null hypothesis that the two 
   readers are interchangeable vs the alternative that there are systematic differences between readers (if chi=TRUE).
   The tests are tests of symmetry (Evans and Hoenig, 1998). If cont.cor=T, then correction for continuity is 
   applied to the McNemar-like chi-square test statistic; the default is to apply the Yates correction but if 
   correct="Edwards" is specified then the correction for continuity is 1.0 instead of 0.5. 
}
\value{Separate lists with tables of various statistics associated with the method.  
}
\references{Evans, G.T. and J.M. Hoenig.  1998. Viewing and Testing Symmetry in Contingency Tables, with 
   Application to Fish Ages. Biometrics 54:620-629.). 
}
\author{John Hoenig, Virginia Institute of Marine Science, 18 December 2012. \email{hoenig@vims.edu}}
\examples{
 data(sbotos)
 compare2(readings=sbotos,usecols=c(1,2),twovsone=TRUE,plot.summary=TRUE,
 barplot=FALSE,chi=TRUE,pool.criterion=1,cont.cor=TRUE,correct="Yates",
 first.name="Reader A",second.name="Reader B")
}
\keyword{misc}