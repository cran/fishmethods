\name{bt.log}
\alias{bt.log}
\title{
Back-transformation of log-transformed mean and variance
}
\description{
Converts a log-mean and log-variance to the original scale and calculates confidence intervals
}
\usage{
bt.log(meanlog = NULL, sdlog = NULL, n = NULL, alpha = 0.05)
}
\arguments{
  \item{meanlog}{sample mean of natural log-transformed values}
  \item{sdlog}{sample standard deviation of natural log-transformed values}
  \item{n}{sample size}
  \item{alpha}{alpha-level used to estimate confidence intervals}
}
\details{
 There are two methods of calcuating the bias-corrected mean on the original scale.
 The \code{bt.mean} is calculated following equation 14 (the infinite series estimation)
 in Finney (1941).  \code{approx.bt.mean} is calculated using the commonly known approximation
 from Finney (1941): 

\emph{mean=exp(meanlog+sdlog^2/2)}.  The variance is 
  \emph{var=exp(2*meanlog)*(Gn(2*sdlog^2)-Gn((n-2)/(n-1)*sdlog^2)} 
  and standard deviation is \emph{sqrt(var)}
 where Gn is  the infinite series function (equation 10).  The variance and standard deviation of the  
  back-transformed mean are \emph{var.mean=var/n; sd.mean=sqrt(var.mean)}. 
 The median is calculated as \emph{exp(meanlog)}.
 Confidence intervals for the back-transformed mean are from the Cox method (Zhou and Gao, 1997) modified
 by substituting the z distribution with the t distribution as recommended by Olsson (2005):

 \emph{LCI=exp(meanlog+sdlog^2/2-t(df,1-alpha/2)*sqrt((sdlog^2/n)+(sdlog^4/(2*(n-1))))} and

 \emph{UCI=exp(meanlog+sdlog^2/2+t(df,1-alpha/2)*sqrt((sdlog^2/n)+(sdlog^4/(2*(n-1))))}

 where \emph{df=n-1}.
}
\value{
A vector containing \code{bt.mean}, \code{approx.bt.mean},\code{var}, \code{sd}, \code{var.mean},\code{sd.mean},
\code{median}, LCI (lower confidence interval), and UCI (upper confidence interval).
}
\references{

Finney, D. J. 1941. On the distribution of a variate whose logarithm is normally distributed. Journal
of  the Royal Statistical Society Supplement 7: 155-161.

Zhou, X-H., and Gao, S. 1997. Confidence intervals for the log-normal mean. Statistics in Medicine 
 16:783-790. 

Olsson, F. 2005. Confidence intervals for the mean of a log-normal distribution. Journal of Statistics 
 Education 13(1). www.amstat.org/publications/jse/v13n1/olsson.html
}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\examples{
## The example below shows accuracy of the back-transformation
y<-rlnorm(100,meanlog=0.7,sdlog=0.2)
known<-unlist(list(known.mean=mean(y),var=var(y),sd=sd(y),
  var.mean=var(y)/length(y),sd.mean=sqrt(var(y)/length(y))))
est<-bt.log(meanlog=mean(log(y)),sdlog=sd(log(y)),n=length(y))[c(1,3,4,5,6)]
known;est}
\keyword{misc}
