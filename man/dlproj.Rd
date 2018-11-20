\name{dlproj}
\alias{dlproj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 This function performs projections for dbsra and catchmsy objects
}
\description{
Make biomass projections by using inputted catch and results of dbsra or catchmsy functions}
\usage{
dlproj(dlobj = NULL, projyears = NULL, projtype = 1, projcatch = NULL, 
grout = 1, grargs = list(lwd = 1, unit = "MT", mains = " ", cex.main = 1, 
cex.axis = 1, cex.lab = 1), grtif = list(zoom = 4, width = 11, height = 13, 
pointsize = 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dlobj}{function dbsra or catchmsy output object}
  \item{projyears}{the number of years for projection.  The first year will be the last year of
	catch data plus one in the original dbsra or catchmsy run.
  }
  \item{projtype}{the type of catch input. 0 = use median MSY from dbsra or catchmsy object, 1 = use
	mean MSY from dbsra or catchmsy object, 2 = user-inputted catch}
  \item{projcatch}{if projtype = 2, a single catch value used over all projection years or a
	vector of catch values (length is equal to \code{projyears}).}
 \item{grout}{numeric argument specifying whether projection graph should be shown on 
  the console only (grout=1) or shown on the console and exported to a TIF graph file (grout=2).
   No graph (grout== 0). If plotted, the median (solid line), mean (dashed line), and 2.5th and 97.5 
   percentiles(dotted lines) are displayed. Use \code{setwd} before running function to direct .tif file
	 to a specific directory. The name of .tif file is automatically determined.} 

 \item{grargs}{list control arguments for plotting functions. \code{lwd} is the line width, \code{unit}
  is the biomass unit for the y-axis label,\code{mains} and 
 \code{cex.main} are the title and character expansion value for the graph, \code{cex.axis} is the 
 character expansion value for the x and y-axis tick labels and \code{cex.lab} is the character 
 expansion value(s) for the x and y-axis labels. 
}
 \item{grtif}{list control arguments for the .TIF graph file. See \code{tiff} help file in R.} 
}
\details{
The biomass estimate of the last year+1 is used as the starting biomass (year 1 in projections)
and leading parameters from each plausible (accepted) run are used to project biomass ahead \code{projyears} years
 using either the MSY estimate (median or mean) from all plausible runs or inputted catch values.
The biomass estimates are loaded from either the "Biotraj-dbsra.csv" or "Biotroj-cmsy.csv" files that were 
automatically saved in functions "dbsra" and "catchmsy".

Use \code{setwd()} before running the function to change the directory where .csv files are stored.  

  
}
\value{
\item{type}{object projection type}
\item{ProjBio}{dataframe of biomass projections for each plausible run}
}
\references{

Martell, S. and R. Froese. 2012. A simple method for estimating MSY from catch and resilience.
Fish and Fisheries 14:504-514.

Dick, E. J. and A. D. MacCall. 2011. Depletion-based stock reduction analysis: a catch-based method for determining
sustainable yield for data-poor fish stocks. Fisheries Research 110: 331-341.
}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\seealso{\code{\link{catchmsy}} \code{\link{dbsra}}}
\examples{
\dontrun{
  data(lingcod)
   outs<-catchmsy(year=lingcod$year,
    catch=lingcod$catch,catchCV=NULL,
    catargs=list(dist="none",low=0,up=Inf,unit="MT"),
    l0=list(low=0.8,up=0.8,step=0),
    lt=list(low=0.01,up=0.25,refyr=2002),sigv=0,
    k=list(dist="unif",low=4333,up=433300,mean=0,sd=0),
    r=list(dist="unif",low=0.015,up=0.5,mean=0,sd=0),
    bk=list(dist="unif",low=0.5,up=0.5,mean=0,sd=0),
    M=list(dist="unif",low=0.24,up=0.24,mean=0.00,sd=0.00),
    nsims=30000)
   outbio<-dlproj(dlobj = outs, projyears = 20, projtype = 0, grout = 1)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
