\name{gap}
\alias{gap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tukey's Gapping
}
\description{
This function finds unusual spaces or gaps in a vector of random samples}
\usage{
gap(x = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  vector of values
  }
 }
\details{
Values (x) are sorted from smallest to largest. Then Z values are calculated as follows:
Zn-i+1=[i*(n-i)(Xn-i+1 - Xn-i)]^0.5

where n is the sample size 

for i = 2,...,n calulate the 25 percent trimmed mean and divide into Z.
This standardizes the distribution of the weighted gaps around a middle value of one. Suspiciously large
observations should correspond to large standardized weighted gaps.
 
}
\value{vector of standardized weighted gaps }
\references{
Tukey, J. W. 1971. Exploratory data analysis. Addison-Wesley, Reading, MA. 431 pp.}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
 y<-c(rnorm(10,10,2),1000)
 gap(y)
}
\keyword{misc}
