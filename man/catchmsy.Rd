\name{catchmsy}
\alias{catchmsy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimating MSY from catch and resilience
}
\description{
This function estimates MSY following Martell and Francis (2012) using the
 generalized Pella-Tomlinson production equation.
}
\usage{
catchmsy(year = NULL, catch = NULL, catchCV = NULL, 
catargs = list(dist = "none", low = 0, up = Inf, unit = "MT"), 
l0 = list(low = 0, up = 1, step = 0), lt = list(low = 0, up = 1, termyr = NULL), 
sigv = 0, k = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), 
r = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), 
bk = list(dist = "unif", low = 0.5, up = 0.5, mean = 0, sd = 0), 
M = list(dist = "unif", low = 0.2, up = 0.2, mean = 0, sd = 0), 
nsims = 10000, grout = 1, graphs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 
grargs = list(lwd = 1, pch = 16, cex = 1, nclasses = 20, mains = " ", 
cex.main = 1, cex.axis = 1, cex.lab = 1), pstats = list(ol = 1, mlty = 1, 
mlwd = 1.5, llty = 3, llwd = 1, ulty = 3, ulwd = 1), grtif = list(zoom = 4, 
width = 11, height = 13, pointsize = 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{vector containing the time series of numeric year labels.}
  \item{catch}{vector containing the time series of catch data (in weight). Missing values are not
   allowed.}
  \item{catchCV}{vector containing the time series of coefficients of variation associated with catch 
           if resampling of catch is desired; otherwise, catchCV = NULL. }
  \item{catargs}{list arguments associated with resampling of catch. \code{dist} is the specification 
      of the resampling distribution to use ("none" = no resampling, "unif"=uniform, "norm" = normal, 
      and "lnorm" =log-normal). If "lnorm" is selected, catch is log transformed and standard deviation
      on the log scale is calculated from the specificed CVs using the relationship sdlog=sqrt(log(CV^2+1)). 
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
  \item{lt}{list arguments for the depletion level in the selected terminal year (\code{termyr}). 
   \code{low} and \code{up} are the lower and upper bounds of depletion level in \code{termyear}. 
    \code{termyear} can range from the first year to the year after the last year of catch (t+1). 
  }
  \item{sigv}{standard deviation of the log-normal random process error.  \code{signv} = 0 for no 
   process error.
  }
  \item{k}{list arguments for the carrying capacity. \code{dist} is the statistical distribution name
   from which to sample \code{k}. \code{low} and \code{up} are the lower and upper bounds of \code{k} 
   in the selected distribution. 
   \code{mean} and \code{sd} are the mean and standard deviation for selected distributions. The 
    following are valid distributions: "unif" - uniform, "norm" - normal, "lnorm" - log-normal, 
    "gamma" - gamma, and "beta" - beta distributions. "unif" requires non-missing values for \code{low} 
    and \code{up}. "norm", "lnorm", "gamma" and "beta", require non-missing values for 
    \code{low},\code{up}, \code{mean} and \code{sd}. If "lnorm" is used, \code{mean} and \code{sd} 
    must be on the natural log scale (keep \code{low} and \code{up} on the original scale). If use 
    of a constant value is desired, specify the "unif" distribution and enter the same value for 
    \code{low} and \code{up} arguments.
  }
\item{r}{list arguments for the intrinsic growth rate. \code{dist} is the statistical distribution name 
  from which to sample \code{r}. \code{low} and \code{up} are the lower and upper bounds of \code{r} 
  in the selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for 
  selected distributions. Valid distributions are the same as in \code{k}.If use of a constant value is
     desired, specify the "unif" distribution and enter the same value for the \code{low} and \code{up}
     arguments.
  }
\item{bk}{list arguments for Bmsy/k. \code{dist} is the statistical distribution name from which 
   to sample \code{Bmsy/k}. \code{low} and \code{up} are the lower and upper bounds of \code{Bmsy/k} in
  the selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for selected 
  distributions. Valid distributions are the same as in \code{k}. Bmsy/k is used in the Pella-Tomlinsom 
  production function to generate the shape parameter (p). If use of a constant value is desired, 
  specify the "unif" distribution and enter the same value for the \code{low} and \code{up} arguments. 
  For the Schaefer production curve, use Bmsy/k = 0.5. 
  }
\item{M}{list arguments for natural mortality. \code{dist} is the statistical distribution name from 
 which to sample \code{M}. \code{low} and \code{up} are the lower and upper bounds of \code{M} in the 
 selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for selected 
 distributions. Valid distributions are the same as in \code{k}. M is used to determine exploitation 
 rate (Umsy) at MSY. If use of a constant value is desired, specify the "unif" distribution and enter 
 the same value for the \code{low} and \code{up} arguments.
  }
\item{nsims}{number of Monte Carlos samples.} 
\item{grout}{numeric argument specifying whether graphs should be printed to console only (1) or to 
 both the console and TIF graph files (2).} 
\item{graphs}{vector specifying which graphs should be produced. 1 = line plot of observed catch versus
  year,2 = point plot of plausible \code{k} versus \code{r} values, 3 = histogram of plausible r values, 
 4 = histogram of plausible k values, 5 = histogram of plausible Bmsy/k values, 6 = histogram of M values,
 7 = histogram of MSY from plausible values of l0,k,r, and Bmsy/k, 8 = histogram of Bmsy from plausible 
 values of l0,k,r, and Bmsy/k, 9 = histogram of Fmsy from plausible values of l0,k,r, and Bmsy/k, 10 = 
histogram of Umsy values from Fmsy and M, 11 = histogram of overfishing limit (OFL) values from Umsys. 
Any combinations of graphs can be selected within c().  Default is all.
  }
\item{grargs}{list control arguments for plotting functions. \code{lwd} is the line width for graph = 1, 
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
information on resilience is known. In this implementation, the Pella-Tomlinson production equation is
 used and an additional parameter (Bmsy/k) can be explored along with to r and k.

The  Pella-Tomlinson production model (Haddon, 2001) is 
 
    B[t+1]<-B[t]+(r/p)*B[t]*(1-(B[t]/k)^p)-catch[t]

where B is biomass in year t, \code{r} is the instrince rate of increase, p is a shape parameter 
(p=1 for Schaefer production model), \code{k} is the carrying capacity and catch is the catch in year 
t. Biomass is the first year is calculated by B[1]=\code{k}*\code{l0}. For sigv>0, the production 
equation is multiplied by exp(rnorm(1,0,sigv)) if process error is desired.  p is calculated by 
specifying the Bmsy/k ratio and iteratively solving for p in the equation:

 Bmsy/k=(1/(1+p))^(1/p)

The maximum sustainable yield (MSY) is calculated as

  MSY=r*k/((1+p)^((p+1)/p))

Biomass at MSY is calculated as

  Bmsy=k*/((1+p)^(1/p))

Fishing mortality at MSY is calculated as 

 Fmsy=r*(1/(1+p))

Exploitation rate at MSY is calculated as

 Umsy=(Fmsy/(Fmsy+M))*(1-exp(-Fmsy-M))

The overfishing limit in \code{termyr} is calculated as

  OFL=B[termyr]*Umsy
}
\value{
\item{Initial}{dataframe containing the initial values for each explored parameter.}
  \item{Parameters}{dataframe containing the mean, median, 2.5th and 97.5% percentiles of the 
        plausible (likelihood=1) parameters. }  
  \item{Estimates}{dataframe containing the mean, median, 2.5th and 97.5% percentiles 
        of the management quantities (i.e., MSY, Bmsy, etc.) from the plausible parameters 
        (likelihood=1)}
  \item{Data}{dataframe containing the values of l0, k, r, Bmsy/k, M and associated management 
       quantities for all (likelihood=0 and likelihood=1) random draws.}
}
\references{

Haddon, M. 2001. Modelling and quantitative methods in fisheries. Chapman & Hall/CRC.

Martell, S. and R. Froese. 2012. A simple method for estimating MSY from catch and resilience.
Fish and Fisheries 14:504-514.
}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}
}
\note{The random distribution function was adapted from Nadarajah, S and S. Kotz. 2006. R 
 programs for computing truncated distributions. Journal of Statistical Software 16, code snippet 2.
}
\examples{
 \dontrun{
   data(lingcod)
   catchmsy(year=lingcod$year,
    catch=lingcod$catch,catchCV=NULL,
    catargs=list(dist="none",low=0,up=Inf,unit="MT"),
    l0=list(low=0.8,up=0.8,step=0),
    lt=list(low=0.01,up=0.25,termyr=2002),sigv=0,
    k=list(dist="unif",low=4333,up=433300,mean=0,sd=0),
    r=list(dist="unif",low=0.015,up=0.5,mean=0,sd=0),
    bk=list(dist="unif",low=0.5,up=0.5,mean=0,sd=0),
    M=list(dist="unif",low=0.24,up=0.24,mean=0.00,sd=0.00),
    nsims=10000)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
