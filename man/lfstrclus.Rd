\name{lfstrclus}
\alias{lfstrclus}
\title{Statistical Comparison of Length Frequencies from Stratified Random Cluster Sampling}
\description{Statistical comparison of length frequencies is performed using the two-sample
 Kolmogorov & Smirnov test. Randomization procedures are used to derive the null probability
 distribution. 
}
\usage{
lfstrclus(group = NULL, strata = NULL, weights = NULL,
 haul = NULL, len = NULL, number = NULL, binsize = NULL,
 resamples = 100)
}
\arguments{
 \item{group}{ vector containing the identifier used for group membership of length data.  This variable is
used to determine the number of groups and comparisons.  Identifier can be numeric or character.}
 \item{strata}{vector containing the numeric identifier used for strata membership of length data.  
   There must be a unique identifier for each stratum regardless of group membership.}
 \item{weights}{vector containing the strata weights (e.g., area, size, etc.) used to calculate the stratified mean length for a group.  
  }
 \item{haul}{vector containing the variable used to identify the sampling unit (e.g., haul) of length data. Identifier can be numeric or character.}
 \item{len}{vector containing the length class. Each length class record must have associated group, strata, weights, and haul identifiers. }
 \item{number}{vector containing the number of fish in each length class.}
 \item{binsize}{size of the length class (e.g., 5-cm, 10, cm, etc.) used to construct the cumulative length frequency
from raw length data.  The formula used to create bins is 
          \eqn{trunc(len/binsize)*binsize+binsize/2}. If use of the raw length classes is desired, then \code{binsize=0}.}
 \item{resamples}{number of randomizations. Default = 100.}
}
\details{
Length frequency distributions of fishes are commonly tested for differences among groups (e.g., regions, sexes, etc.) using a two-sample Kolmogov-Smirnov test (K-S).  Like most statistical
tests, the K-S test requires that observations are collected at random and are independent of each other to satisfy assumptions. These basic assumptions are
violated when gears (e.g., trawls, haul seines, gillnets, etc.) are used to sample fish because individuals are collected in clusters .  In this case, the "haul", not the individual fish, is the
 primary sampling unit and statistical comparisons must take this into account.  

To test for difference between length frequency distributions from stratified random cluster sampling, a randomization test that uses "hauls" as 
the primary sampling unit can be used to generate the null probability distribution.  In a randomization test, an observed test statistic is compared to an empirical 
probability density distribution of a test statistic under the null hypothesis of no difference.  The observed test statistic
used here is the Kolmogorov-Smirnov statistic (Ds) under a two-tailed test:

   					 \deqn{Ds= max|S1(X)-S2(X)|}

where S1(X) and S2(X) are the observed cumulative proportions at length for group 1 and group 2 in the paired comparisons.

Proportion of fish of length class j in strata-set (group variable) used to derive \code{Ds} is calculated as

\deqn{p_j=\frac{\sum{A_k\bar{X}}_{jk}}{\sum{A_k\bar{X}}_k}}{%
      p(j)=sum(Ak Xjk)/sum(Ak Xk)}

where \eqn{A_k} is the weight of stratum k, \eqn{\bar{X}_{jk}} is the mean number per haul of length class \code{j} in stratum \code{k}, and
\eqn{\bar{X}_k}  is the mean number per haul in stratum \code{k}. The numerator and denominator are summed over all \code{k}. Before calculation of
cumulative proportions, the length class distributions for each group are corrected for missing lengths and are
constructed so that the range and intervals of each distribution match. 

It is assumed all fish caught are measured. If subsampling occurs, the numbers at length (measured) must be expanded to the total caught.
 
To generate the empirical probability density function (pdf), length data of hauls from all strata are pooled and then hauls are randomly assigned without replacement
to each stratum with haul sizes equal to the original number of stratum hauls.  Cumulative proportions are
then calculated as described above. The K-S statistic is calculated from the cumulative length frequency distributions of the two groups 
of randomized data.  The randomization procedure is repeated \code{resamples} times to 
obtain the pdf of D.  To estimate the significance of Ds, the proportion of all randomized D values
that were greater than or equal to Ds, including the occurrence of Ds in the numerator and denominator, is calculated
(Manly, 1997).
  
Data vectors described in \code{arguments} should be aggregated so that each record contains the number of fish in each length class by group, strata, weights, and haul identifier. For example,


\tabular{lrrrrr}{
\code{group} \tab \code{stratum} \tab \code{weights} \tab \code{tow} \tab \code{length} \tab \code{number} \cr
North	\tab 10 \tab 88 \tab 1 \tab 10 \tab 2 \cr
North	\tab 10 \tab 88 \tab 1 \tab 12 \tab 5 \cr
North	\tab 10 \tab 88 \tab 2 \tab 11 \tab 3 \cr
North	\tab 11 \tab 103 \tab 1 \tab 10 \tab 17 \cr
North	\tab 11 \tab 103 \tab 2 \tab 14 \tab 21 \cr
.      \tab . \tab .  \tab . \tab .  \tab . \cr
.      \tab . \tab .  \tab . \tab .  \tab . \cr
South	\tab 31 \tab 43 \tab 1 \tab 12 \tab 34 \cr
South	\tab 31 \tab 43 \tab 1 \tab 14 \tab 3
}

To correctly calculate the stratified mean number per haul, zero tows must be included in the dataset. 
To designate records for zero tows, fill the length class and number at length with zeros.  The first line in
the following table shows the appropriate coding for zero tows:

\tabular{lrrrrr}{
\code{group} \tab \code{stratum} \tab \code{weights} \tab \code{tow} \tab \code{length} \tab \code{number} \cr
North	\tab 10 \tab 88 \tab 1 \tab 0 \tab 0 \cr
North	\tab 10 \tab 88 \tab 2 \tab 11 \tab 3 \cr
North	\tab 11 \tab 103 \tab 1 \tab 10 \tab 17 \cr
North	\tab 11 \tab 103 \tab 2 \tab 14 \tab 21 \cr
.      \tab . \tab .  \tab . \tab .  \tab . \cr
.      \tab . \tab .  \tab . \tab .  \tab . \cr
South	\tab 31 \tab 43 \tab 1 \tab 12 \tab 34 \cr
South	\tab 31 \tab 43 \tab 1 \tab 14 \tab 3
}
		

}
\value{
  \item{results}{list element containing the Ds statistics from the observed data comparisons and significance probabilities.}
 \item{obs_prop}{list element containing the cumulative proportions from each group.} 
 \item{Drandom}{list element containing the D statistics from randomization for each comparison.}
}
\references{ Manly, B. F. J. 1997.  Randomization, Bootstrap and Monte Carlos Methods in Biology.
Chapman and Hall, New York, NY, 399 pp.

Seigel, S. 1956. Nonparametric Statistics for Behavioral Sciences. McGraw-Hill, New York, NY. 312 p.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{lfclus}}}
\examples{
data(codstrcluslen)
lfstrclus(
group=codstrcluslen$region,strata=codstrcluslen$stratum,
 weights=codstrcluslen$weights,haul=codstrcluslen$tow,
 len=codstrcluslen$length,number=codstrcluslen$number,
 binsize=5,resamples=100)
}
\keyword{htest}
\keyword{misc}

