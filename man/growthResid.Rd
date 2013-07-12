\name{growthResid}
\alias{growthResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot residuals of growth model fitted to tag data}
\description{
Plot residuals (observed - expected growth increments) vs relative
  age at the time of tagging and versus time at liberty.}
\usage{
growthResid(K, Linf, dat, lentag, lenrec, timelib, graph =1, 
           main = "Residuals of growth increments",
           cex.lab=1.5, cex.axis=1.5, cex.main=1,
           xlab1="Relative age, yr", xlab2="Time at liberty, yr",
           ylab="Observed - expected increment",
           xlim1=NULL, xlim2=NULL, ylim=NULL, col=1, returnvec=FALSE, 
           returnlimits=FALSE, warn=TRUE,...) }
\arguments{
  \item{K}{parameter of the von Bertalanffy growth equation}
  \item{Linf}{parameter of the von Bertalanffy growth equation}
  \item{dat}{dataframe containing length at tagging, length at recapture
                   and time at liberty. These must be named lentag, lenrec
                   and timelib or else column 1 must contain the length at
                   tagging, column 2 must contain length at recapture and
                   column 3 must contain time at liberty}
  \item{lentag}{alternative way to pass data to function}
  \item{lenrec}{alternative way to pass data to function}
  \item{timelib}{alternative way to pass data to function}
  \item{graph}{which graph to plot - 1: residuals versus Relative age, 2: residuals versus time-at-liberty}
  \item{main}{an overall title for the plot}
  \item{cex.lab}{The magnification to be used for x and y labels relative to
          the current setting of cex}
  \item{cex.axis}{The magnification to be used for axis annotation relative to 
           the current setting of cex}
 \item{cex.main}{The magnification to be used for main titles relative to 
            the current setting of cex}
 \item{xlab1}{a title for the x axis 1}
 \item{xlab2}{a title for the x axis 2}
 \item{ylab}{a title for the y axis}
 \item{xlim1}{lower and upper limits of x axis 1 e.g., c(0,100)}
 \item{xlim2}{lower and upper limits of x axis 2 e.g., c(0,100)}
 \item{ylim}{lower and upper limits of y axis e.g., c(0,100)}
 \item{col}{color of points in plot}
 \item{returnvec}{logical - if TRUE, function returns a dataframe with the
        computed age at tagging and the residual (obs - pred increment)}
 \item{returnlimits}{logical - if TRUE, function returns the x and y limits for the plot}
 \item{warn}{logical - if TRUE, function issues a warning if names of variables
                   in dat do not match the 3 names expected.}
 \item{\dots}{other arguments to pass to \emph{plot}}

}
\details{
  Function plots residuals (observed - expected growth increments) vs relative
  age at the time of tagging and vs time at liberty from VB growth model fitted 
  to tagging data. Relative age is calculated by inverting the von Bertalanffy growth curve. 
}
\value{
  \item{output}{If returnvec = TRUE, computed age and residuals. If returnlimits=TRUE, 
       x and y limits for plot}
}
\author{
Janos Hoenig Virginia Institute of Marine Science May 2013 \email{hoenig@vims.edu }
}
\examples{
 data(bonito)
 temp<-bonito[c(bonito$T2-bonito$T1)>0,]
 growthResid(0.19,97.5,lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2-temp$T1),graph=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
