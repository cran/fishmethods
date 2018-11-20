\name{clus.lf}
\alias{clus.lf}
\title{Statistical Comparison of Length Frequencies from Simple Random Cluster Sampling}
\description{Statistical comparison of length frequencies is performed using the two-sample
 Kolmogorov & Smirnov test. Randomization procedures are used to derive the null probability
 distribution. 
}
\usage{
clus.lf(group = NULL, haul = NULL, len = NULL, number= NULL,
 binsize = NULL, resamples = 100)
}
\arguments{
  \item{group}{ vector containing the identifier used for group membership of length data.  This variable is
used to determine the number of groups and comparisons.  Identifier can be numeric or character.}
 \item{haul}{vector containing the variable used to identify the sampling unit (e.g., haul) of length data. Identifier can be numeric or character.}
  \item{len}{vector containing the length class data. There should be one record for each length class by group and haul. }
 \item{number}{vector containing the numbers of fish in each length class.}
  \item{binsize}{size of the length class (e.g., 5-cm, 10, cm, etc.) used to construct the cumulative length frequency
from raw length data.  The formula used to create bins is \eqn{trunc(len/binsize)*binsize+binsize/2}. 
If use of the raw length classes is desired, then \code{binsize=0}.}
  \item{resamples}{number of randomizations. Default = 100.}
}

\details{
Length frequency distributions of fishes are commonly tested for differences among groups (e.g., regions, sexes, etc.) using a two-sample Kolmogov-Smirnov test (K-S).  Like most statistical
tests, the K-S test requires that observations are collected at random and are independent of each other to satisfy assumptions. These basic assumptions are
violated when gears (e.g., trawls, haul seines, gillnets, etc.) are used to sample fish because individuals are collected in clusters .  In this case, the "haul", not the individual fish, is the
 primary sampling unit and statistical comparisons must take this into account.  

To test for difference between length frequency distributions from simple random cluster sampling, a randomization test that uses "hauls" as 
the primary sampling unit can be used to generate the null probability distribution.  In a randomization test, an observed test statistic is compared to an empirical 
probability density distribution of a test statistic under the null hypothesis of no difference.  The observed test statistic
used here is the Kolmogorov-Smirnov statistic (Ds) under a two-tailed test:


   					 \deqn{Ds= max|S1(X)-S2(X)|}

where S1(X) and S2(X) are the observed cumulative length frequency distributions of group 1 and group 2 in the paired comparisons. 
S1(X) and S2(X) are calculated such that \code{S(X)=K/n} where K is the number of scores equal to or less
than X and n is the total number of length observations (Seigel, 1956).  

To generate the empirical probability density function (pdf), haul data are randomly assigned without replacement to the two groups with samples sizes equal to the original number of hauls in each group under comparison.
The K-S statistic is calculated from the cumulative length frequency distributions of the two groups 
of randomized data.  The randomization procedure is repeated \code{resamples} times to 
obtain the pdf of D.  To estimate the significance of Ds, the proportion of all randomized D values
that were greater than or equal to Ds is calculated.
  
It is assumed all fish caught are measured. If subsampling occurs, the number at length (measured) must be expanded to the total caught.
 
Data vectors described in \code{arguments} should be aggregated so that each record contains the number of fish in each length class by group and haul identifier. For example,

\tabular{rrrr}{
\code{group} \tab \code{tow} \tab \code{length} \tab \code{number} \cr
North	\tab 1 \tab 10 \tab 2 \cr
North	\tab 1 \tab 12 \tab 5 \cr
North	\tab 2 \tab 11 \tab 3 \cr
North	\tab 1 \tab 10 \tab 17 \cr
North	\tab 2 \tab 14 \tab 21 \cr
.     \tab . \tab .  \tab . \cr
.     \tab . \tab .  \tab . \cr
South	\tab 1 \tab 12 \tab 34 \cr
South	\tab 1 \tab 14 \tab 3
}

}
\value{
\item{results}{list element containing the Ds statistics from the observed data comparisons and significance probabilities.}
\item{obs_prop}{list element containing the observed cumulative proportions for each group.}
  \item{Drandom}{list element containing the D statistics from randomization for each comparison.}
}
\references{ Manly, B. F. J. 1997.  Randomization, Bootstrap and Monte Carlos Methods in Biology.
Chapman and Hall, New York, NY, 399 pp.

Seigel, S. 1956. Nonparametric Statistics for Behavioral Sciences. McGraw-Hill, New York, NY. 312 p.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{clus.str.lf}}}
\examples{
data(codcluslen)
clus.lf(group=codcluslen$region,haul=c(paste("ST-",codcluslen$tow,sep="")),
 len=codcluslen$length, number=codcluslen$number, 
 binsize=5,resamples=100)
}
\keyword{misc}
