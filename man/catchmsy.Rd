\name{catchmsy}
\alias{catchmsy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimating MSY from catch and resilience
}
\description{
This function estimates MSY following Martell and Froese(2012).
}
\usage{
catchmsy(year = NULL, catch = NULL, catchCV = NULL, 
catargs = list(dist = "none", low = 0, up = Inf, unit = "MT"), 
l0 = list(low = 0, up = 1, step = 0), lt = list(low = 0, up = 1, 
refyr = NULL), 
sigv = 0, k = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), 
r = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), 
M = list(dist = "unif", low = 0.2, up = 0.2, mean = 0, sd = 0), 
nsims = 10000, catchout = 0, grout = 1, 
graphs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 
grargs = list(lwd = 1, pch = 16, cex = 1, nclasses = 20, mains = " ", 
cex.main = 1, 
cex.axis = 1, cex.lab = 1), 
pstats = list(ol = 1, mlty = 1, mlwd = 1.5, llty = 3, llwd = 1, ulty = 3, 
ulwd = 1), 
grtif = list(zoom = 4, width = 11, height = 13, pointsize = 10))
}
\arguments{
  \item{year}{vector containing the time series of numeric year labels.}
  \item{catch}{vector containing the time series of catch data (in weight). Missing values are not
   allowed.}
  \item{catchCV}{vector containing the time series of coefficients of variation associated with catch 
           if resampling of catch is desired; otherwise, catchCV = NULL. }
  \item{catargs}{list arguments associated with resampling of catch. \code{dist} is the specification 
      of the resampling distribution to use ("none" = no resampling, "unif"=uniform, "norm" = normal, 
      and "lnorm" =log-normal). If "lnorm" is selected, catch is log transformed and standard deviation
      on the log scale is calculated from the specified CVs using the relationship sdlog=sqrt(log(CV^2+1)). 
       \code{low} and \code{up} are the lower and upper limit of distribution (if truncation is desired). 
       \code{unit} is the weight unit of catch (used in graph labels; default="MT"). If "unif", the
       \code{catch} must be incorporated in \code{low} and \code{up} arguments.  For instance, if the 
       lower limit to sample is the value of \code{catch}, specify \code{low}=catch. If some maximum 
       above catch will be the upper limit, specify \code{up}=50*catch.  The limits for each year will 
       be applied to catch internally. 
     }
  \item{l0}{list arguments for the relative biomass in year 1. \code{low} and \code{up} are the lower 
   and upper bounds of the starting value of relative biomass (in relation to k) in year 1. \code{step} 
   is the step increment to examine.  If \code{step}=0, then \code{l0} is randomly selected from a
   uniform distribution using the lower and upper starting values. If \code{step}>0, then step increments
   are used (in this case, the number of simulations (\code{nsims}) are used for each increment). 
  }
  \item{lt}{list arguments for the depletion level in the selected reference year (\code{refyr}). 
   \code{low} and \code{up} are the lower and upper bounds of depletion level in \code{refyr}. 
    \code{refyr} can range from the first year to the year after the last year of catch (t+1). 
  }
  \item{sigv}{standard deviation of the log-normal random process error.  \code{signv} = 0 for no 
   process error.
  }
  \item{k}{list arguments for the carrying capacity. \code{dist} is the statistical distribution name
   from which to sample \code{k}. \code{low} and \code{up} are the lower and upper bounds of \code{k} 
   in the selected distribution. 
   \code{mean} and \code{sd} are the mean and standard deviation for selected distributions. The 
    following are valid distributions: "none", "unif" - uniform, "norm" - normal, "lnorm" - log-normal, 
    "gamma" - gamma, and "beta" - beta distributions. "unif" requires non-missing values for \code{low} 
    and \code{up}. "norm", "lnorm", "gamma" and "beta", require non-missing values for 
    \code{low},\code{up}, \code{mean} and \code{sd}. If "lnorm" is used, \code{mean} and \code{sd} 
    must be on the natural log scale (keep \code{low} and \code{up} on the original scale). If 
    \code{dist} = "none", the mean is used as a fixed value.
  }
\item{r}{list arguments for the intrinsic growth rate. \code{dist} is the statistical distribution name 
  from which to sample \code{r}. \code{low} and \code{up} are the lower and upper bounds of \code{r} 
  in the selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for 
  selected distributions. Valid distributions are the same as in \code{k}. If \code{dist} = "none", 
  the mean is used as a fixed value.

  }
\item{M}{list arguments for natural mortality. \code{dist} is the statistical distribution name from 
 which to sample \code{M}. \code{low} and \code{up} are the lower and upper bounds of \code{M} in the 
 selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for selected 
 distributions. Valid distributions are the same as in \code{k}. M is used to determine exploitation 
 rate (Umsy) at MSY. If \code{dist} = "none", the mean is used as a fixed value.
  }
\item{nsims}{number of Monte Carlos samples.} 
\item{catchout}{If resampling \code{catch}, save catch trajectories from each Monte Carlos simulation 
 - 0 = No (default), 1 = Yes.} 
\item{grout}{numeric argument specifying whether graphs should be printed to console only (1) or to 
 both the console and TIF graph files (2).Use \code{setwd} before running function to direct .tif files
 to a specific directory. Each name of each file is automatically determined.} 
\item{graphs}{vector specifying which graphs should be produced. 1 = line plot of observed catch versus
  year,2 = point plot of plausible \code{k} versus \code{r} values, 3 = histogram of plausible r values, 
 4 = histogram of plausible k values,  5 = histogram of M values,
 6 = histogram of MSY from plausible values of l0,k,r, and Bmsy/k, 7 = histogram of Bmsy from plausible 
 values of l0,k,r, and Bmsy/k, 8 = histogram of Fmsy from plausible values of l0,k,r, and Bmsy/k, 9 = 
histogram of Umsy values from Fmsy and M, 10 = histogram of overfishing limit (OFL) in last year+1 values 
from Umsys, and 11 = line plots of accepted and rejected biomass trajectores with median and 2.5th and 97.5th
percentiles (in red).  Any combinations of graphs can be selected within c().  Default is all.
  }
\item{grargs}{list control arguments for plotting functions. \code{lwd} is the line width for graph = 1 and 11, 
 \code{pch} and \code{cex} are the symbol character and character expansion value used in graph = 2, 
 \code{nclasses} is the nclass argument for the histogram plots (graphs 3-11), \code{mains} and 
 \code{cex.main} are the titles and character expansion values for the graphs, \code{cex.axis} is the 
 character expansion value(s) for the x and y-axis tick labels and \code{cex.lab} is the character 
 expansion value(s) for the x and y-axis labels.  Single values of \code{nclasses},\code{mains}, 
 \code{cex.main},\code{cex.axis}, \code{cex.lab} are applied to all graphs.  To change arguments for 
 specific graphs, enclose arguments within c() in order of the number specified in \code{graphs}. 
}
\item{pstats}{list control arguments for plotting the mean and 95% confidence intervals of parameters 
 and management quantities on respective graphs. \code{ol} = 0, do not overlay values on plots, 1 = 
 overlay values on plots. \code{mlty} and \code{mlwd} are the line type and line width of the mean value;
  \code{llty} and \code{llwd} are the line type and line wdith of the 2.5% percentile; code{ulty} and 
\code{ulwd} are the line type and line width of the 97.5% percentile.
}
\item{grtif}{list arguments for the .TIF graph files. See \code{tiff} help file in R.} 
}
\details{
The method of Martell and Froese (2012) is used to produce estimates of MSY where only catch and 
information on resilience is known. 

The Schaefer production model is 
 
    B[t+1]<-B[t]+r*B[t]*(1-B[t]/k)-catch[t]

where B is biomass in year t, \code{r} is the instrince rate of increase, \code{k} is the carrying 
capacity and \code{catch} is the catch in year t. Biomass is the first year is calculated by 
B[1]=\code{k}*\code{l0}. For sigv>0, the production equation is multiplied by exp(rnorm(1,0,sigv)) 
if process error is desired. 
The maximum sustainable yield (MSY) is calculated as

  MSY=r*k/4

Biomass at MSY is calculated as

  Bmsy=k/2

Fishing mortality at MSY is calculated as 

 Fmsy=r/2

Exploitation rate at MSY is calculated as

 Umsy=(Fmsy/(Fmsy+M))*(1-exp(-Fmsy-M))

The overfishing limit in last year+1  is calculated as

  OFL=B[last year +1]*Umsy

\code{length(year)+1} biomass estimates are made for each run.

If using the R Gui (not Rstudio), run 

      graphics.off()
  	windows(width=10, height=12,record=TRUE)
      .SavedPlots <- NULL

before running the catchmsy function to recall plots.
}
\value{
  \item{Initial}{dataframe containing the initial values for each explored parameter.}
  \item{Parameters}{dataframe containing the mean, median, 2.5th and 97.5% percentiles for the 
        plausible (likelihood=1) parameters. }  
  \item{Estimates}{dataframe containing the mean, median, 2.5th and 97.5% percentiles 
        of the management quantities (i.e., MSY, Bmsy, etc.) for the plausible parameters 
        (likelihood=1)}
  \item{Values}{dataframe containing the values of l0, k, r, Bmsy/k, M and associated management 
       quantities for all (likelihood=0 and likelihood=1) random draws.}
  
 \item{end1yr}{value of the last year of catch data + 1 for use in function \code{dlproj}.}
 \item{type}{designates the output object as a \code{catchmsy} object for use in function \code{dlproj}.}  

 The biomass estimates from each simulation are not stored in memory but are automatically 
     saved to a .csv file named "Biotraj-cmsy.csv". Yearly values for each simulation are stored across 
     columns.  The first column holds the likelihood values for each simulation (1= accepted, 0 = rejected).  
     The number of rows equals the number of simulations (\code{nsims}). This file is loaded to plot
      graph 11 and it must be present in the default or \code{setwd()} directory.  
    
  When \code{catchout}=1,   catch values randomly selected are saved to a .csv file named "Catchtraj-cmsy.csv". 
     Yearly values for each simulation are stored across columns.  The first column holds the likelihood 
     values (1= accepted, 0 = rejected).  The number of rows equals the number of simulations (\code{nsims}).   
  
 Use \code{setwd()} before running the function to change the directory where .csv files are stored.  
}
\references{

Martell, S. and R. Froese. 2012. A simple method for estimating MSY from catch and resilience.
Fish and Fisheries 14:504-514.
}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}
}
\note{The random distribution function was adapted from Nadarajah, S and S. Kotz. 2006. 
 R programs for computing truncated distributions. Journal of Statistical Software 16, 
 code snippet 2.
}
\seealso{\code{\link{dbsra}} \code{\link{dlproj}}}
\examples{
  \dontrun{
   data(lingcod)
   outpt<-catchmsy(year=lingcod$year,
     catch=lingcod$catch,catchCV=NULL,
     catargs=list(dist="none",low=0,up=Inf,unit="MT"),
    l0=list(low=0.8,up=0.8,step=0),
    lt=list(low=0.01,up=0.25,refyr=2002),sigv=0,
    k=list(dist="unif",low=4333,up=433300,mean=0,sd=0),
    r=list(dist="unif",low=0.015,up=0.1,mean=0,sd=0),
    M=list(dist="unif",low=0.18,up=0.18,mean=0.00,sd=0.00),
    nsims=30000)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
