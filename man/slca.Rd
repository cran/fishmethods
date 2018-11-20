\name{slca}
\alias{slca}
\title{A Weakly Parametric Method for the Analysis of Length Composition Data}
\description{
 Shepherd's method for the decomposition of seasonal length frequencies into age classes.
}
\usage{
slca(x, type = 1, fryr=NULL, Linf = NULL, K = NULL, t0 = NULL,
 Lrange = NULL, Krange = NULL)
}
\arguments{
  \item{x}{the dataframe containing the seasonal length frequencies.  The first column contains the lower limit of the
   length bin as a single numeric value, and the second and remaining columns contain the number of fish in each length bin
   for each seasonal length frequency.  The increment of length frequencies should be constant, e.g. every 3 cm.  Empty cells
   must be coded as zeros.  Column headers are not required.}
  \item{type}{the analysis to be conducted: 1=\emph{explore}, 2=\emph{evaluate}.}
  \item{fryr}{the fraction of the year corresponding to when each seasonal length frequency was collected.
Enter one numeric value for each length frequency separated by commas within the concatentation function, e.g. c(0.2,0.45).
Values must be entered for type=1 and type=2. }
  \item{Linf}{the von Bertalanffy L-infinity parameter.  If type=2, then value must be entered.}
  \item{K}{ the von Bertalanffy growth parameter.  If type=2, then value must be entered.}
  \item{t0}{the von Bertalanffy t-sub zero parameter.  If type=2, the value must be entered.}
  \item{Lrange}{the L-infinity range (minimum and maximum) and increment to explore. If type=1, then values must by entered. 
   The first position is the minimum value, the second position is the maximum value, and the third position is the
increment. Values should be separated by commas within the concatentation function, e.g. c(100,120,10).}
  \item{Krange}{the K range and increment to explore. If type=1, then values must by entered. 
The first position is the minimum value, the second position is the maximum value, and the third position is the
increment. Values should be separated by commas within the concatentation function, e.g. c(0.1,0.3,0.02).}
}
\details{
  There are two analytical steps.  In the "explore" analysis, a set of von Bertalanffy parameters that 
best describes the growth of the seasonal length groups is selected from a table of goodness-of-fit measures mapped 
over the range of specified K and L-infinity values.  Once the best K and L-infinity parameters are selected, the corresponding 
t0 value is obtained off the second table. In the "evaluate" analysis, the selected parameters are used to 'slice'
the seasonal length frequencies into age classes.
}
\value{
  If type=1, tables of goodness of fit measures versus L-infinity and K parameters, and t0 values versus L-infinity and K parameters.
  If type=2, table of age classes produced from slicing the length frequencies.
}
\references{ Shepherd, J. G. 1987. A weakly parametric method for the analysis of length composition data. In: D. Pauly and
G. Morgan, (eds). The Theory and Application of Length-Based Methods of Stock Assessment. ICLARM Conf. Ser. Manilla.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\note{Shepherd's Fortran code provided in his original working document was translated into R code. 
}
\examples{

#Data are from Shepherd working document - seasonal length frequencies
# for Raja clavata.
data(Shepherd)

#explore
slca(Shepherd,1,fryr=c(0.2,0.45,0.80),Lrange=c(100,150,10),
Krange=c(0.1,0.3,0.02))

#evaluate
slca(Shepherd,2,fryr=c(0.2,0.45,0.80),Linf=120,K=0.2,t0=0.57)

}
\keyword{misc}
