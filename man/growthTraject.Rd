\name{growthTraject}
\alias{growthTraject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plot growth trajectories obtained from tagging data
}
\description{
 Age and length coordinates for the time of tagging and time of recapture are plotted as line 
 segments overlayed on the von Bertalannfy growth curve}
\usage{
  growthTraject(K, Linf, dat, lentag, lenrec, timelib, subsets=NULL,  
               main = "Growth trajectories & fitted curve",
               cex.lab=1.5, cex.axis=1.5, cex.main=1,
               xlab="Relative age, yr", ylab="Length, cm",
               xlim=NULL, ylim=NULL,ltytraject=1, lwdtraject=1,
               coltraject=1, ltyvonB=1, lwdvonB=2, colvonB="red", 
               returnvec=FALSE, returnlimits=FALSE, warn=TRUE, ...) }
\arguments{
  \item{K}{parameter of the von Bertalanffy growth equation}
  \item{Linf}{parameter of the von Bertalanffy growth equation}
  \item{dat}{dataframe containing length at tagging, length at recapture
                   and time at liberty. These must be named lentag, lenrec
                   and timelib or else column 1 must contain the length at
                   tagging, column 2 must contain length at recapture and
                   column 3 must contain time at liberty OR the variables must 
                   be named lentag, lenrec and timelib}
  \item{lentag}{alternative way to pass data to function}
  \item{lenrec}{alternative way to pass data to function}
  \item{timelib}{alternative way to pass data to function}
  \item{subsets}{factor or integer variable specifying subsets of the data 
    to be plotted with separate colors or line types}
  \item{main}{an overall title for the plot}
  \item{cex.lab}{The magnification to be used for x and y labels relative to
          the current setting of cex}
  \item{cex.axis}{The magnification to be used for axis annotation relative to 
           the current setting of cex}
 \item{cex.main}{The magnification to be used for main titles relative to 
            the current setting of cex}
 \item{xlab}{a title for the x axis}
 \item{ylab}{a title for the y axis}
 \item{xlim}{lower and upper limits of x axis e.g., c(0,100)}
 \item{ylim}{lower and upper limits of y axis e.g., c(0,100)}
 \item{ltytraject}{line type for the growth trajectories}
\item{lwdtraject}{line width for the growth trajectories}
\item{coltraject}{line color for the growth trajectories}
\item{ltyvonB}{line type for the fitted von Bertalanffy growth curve}
\item{lwdvonB}{line width for the fitted von Bertalanffy growth curve}
\item{colvonB}{line color for the fitted von B. curve}
 \item{returnvec}{logical for whether the coordinates of the line
    segments should be returned)}
 \item{returnlimits}{logical for whether the x-axis and y-axis limits
                        should be returned}
 \item{warn}{logical - if TRUE, function issues a warning if names of variables
                   in dat do not match the 3 names expected.}
 \item{\dots}{other arguments to pass to \emph{plot}}

}
\details{
  The relative age at tagging is computed from the inverted von Bertalannfy growth equation (i.e., age expressed as a 
  function of length); the age at recapture is taken to be the age at tagging plus the time at liberty. Then 
  the (age, length) coordinates for the time of tagging and time of recapture are plotted as a line segment. 
  Additional parameters control the format of the plot as follows. A call to plot() sets up the axes. 
  Then a call to arrows() draws the line segments. Finally, a call to curve() adds the von Bertalanffy growth curve. 
  Specifying additional graphical parameters is permissable but these will be passed only to plot().
}
\value{
  \item{output}{if returnvec = TRUE, coordinates of the line segments are returned.
    If returnlimits=TRUE, x and y limits for plot are returned}
}
\author{
Janos Hoenig Virginia Institute of Marine Science May 2013 \email{hoenig@vims.edu }
}
\examples{
 data(bonito)
 temp<-bonito[c(bonito$T2-bonito$T1)>0,]
 growthTraject(0.19,97.5,lentag=temp$L1, lenrec=temp$L2,timelib=c(temp$T2-temp$T1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
