\name{powertrend}
\alias{powertrend}
\title{Power Analysis For Detecting Trends}
\description{
  Power analysis for detecting trends in linear regression is implemented following procedures in Gerrodette (1987; 1991).
}
\usage{
powertrend(trend = 1, A1 = NULL, PSE = NULL, pserel = 1,
 maxyrs = 3, pR = 100, step = 5, alpha = 0.05, tail = 2)
}

\arguments{
  \item{trend}{1 = Linear, 2 = Exponential. Default = 1.}
  \item{A1}{the start year abundance.  In actuality, it can be population size, productivity, diversity,
          mortality rate, etc.}
  \item{PSE}{ the proportional standard error (SE(A)/A) = CV in Gerrodette (1987;1991).}
  \item{pserel}{the relationship between abundance and PSE: 1 = 1/sqrt(A1), 2 = constant, 3 = sqrt(A1). Default = 1.}
  \item{maxyrs}{the maximum number of samples or years to project start year abundance. Default = 3.}
  \item{pR}{ the highest positive percent change to investigate. Default = 100.}
  \item{step}{the increment of the range of percent change to investigate. Default = 5.}
  \item{alpha}{the alpha level (Type I error) to use. Default = 0.05. }
  \item{tail}{type of tailed test: 1 = one-tailed, 2= two-tailed. Default = 2.}
}

\details{
 The probability that an upward or downward trend in abundance (power) will be detected is calculated using linear regression
given number of samples (\code{maxyrs}), estimates of sample variability (\code{PSE}) and abundance-PSE relationship (\code{pserel}),
and percent rate of change.  The program calculates power for each \code{step} increment beginning at -100 percent for declining changes 
and ending at \code{pR} percent for increasing changes. See Gerrodette (1987;1991)
for full details. It is assumed that time intervals between samplings is equal.
}
\value{
  Dataframe containing columns of number of samples (\code{years}), trend selected (\code{trend}), the PSE (\code{pse}), 
  alpha level (\code{alpha}), tail of test (\code{tail}), percent change (\code{R}) over \code{maxyrs}, and power (\code{power}).
}
\references{Gerrodette, T. 1987. A power analysis for detecting trends. Ecology. 68(5): 1364-1372.

Gerrodette, T. 1991. Models for power of detecting trends - a reply to Link and Hatfield. Ecology 72(5): 1889-1892.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
powertrend(A1=1000,PSE=0.1)
}
\keyword{misc}