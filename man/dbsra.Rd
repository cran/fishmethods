\name{dbsra}
\alias{dbsra}
\title{Depletion-Based Stock Reduction Analysis 
}
\description{
This function estimates MSY from catch following Dick and MAcCall (2011).
}
\usage{
dbsra(year = NULL, catch = NULL, catchCV = NULL, 
catargs = list(dist = "none", low = 0, up = Inf, unit = "MT"), 
agemat = NULL, maxn=25, k = list(low = 0, up = NULL, tol = 0.01, permax = 1000), 
b1k = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0),
btk = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0, refyr = NULL), 
fmsym = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), 
bmsyk = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), 
M = list(dist = "unif", low = 0, up = 1, mean = 0, sd = 0), nsims = 10000, 
catchout = 0, grout = 1, 
graphs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), 
grargs = list(lwd = 1, cex = 1, nclasses = 20, mains = " ", cex.main = 1, 
cex.axis = 1, 
cex.lab = 1), pstats = list(ol = 1, mlty = 1, mlwd = 1.5, llty = 3, llwd = 1, 
ulty = 3, ulwd = 1), 
grtif = list(zoom = 4, width = 11, height = 13, pointsize = 10))}
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
  \item{agemat}{median age at entry to the reproductive biomass.}
  \item{maxn}{the maximum limit of the Pella-Tomlinson shape parameter that should not be 
      exceeded in the rule for accepting a run.}  
  \item{k}{list arguments for estimation of  \code{k} (carrying capacity). \code{low} and \code{up} are 
   the lower and upper bounds of the minimization routine and \code{tol} is the tolerance level 
   for minimization. A simple equation \code{((btk)-(b[refyr]/k))^2} is used as
   the objective function. R function \code{optimize} is used to find \code{k}. \code{btk} is described 
   below. \code{permax} is the absolute percent difference between the maximum biomass estimate
   and \code{k} that should not be exceeded in the rule for accepting a run (see \code{details}). 
  }
 \item{b1k}{list arguments for B1/K, the relative depletive level in the first year. 
   \code{dist} is the statistical distribution name from which to sample \code{b1k}. 
   \code{low} and \code{up} are the lower and upper bounds of \code{b1k} 
   in the selected distribution.  \code{mean} and \code{sd} are the mean and standard deviation 
    for selected distributions. The following are valid distributions: "none", "unif" - uniform, 
    "norm" - normal, "lnorm" - log-normal, "gamma" - gamma, and "beta" - beta distributions. 
    "unif" requires non-missing values for \code{low} and \code{up}. "norm", "lnorm", 
    "gamma" and "beta" require non-missing values for \code{low},\code{up}, \code{mean} and 
    \code{sd}. If "lnorm" is used, \code{mean} and \code{sd} must be on the natural log scale
    (keep \code{low} and \code{up} on the original scale). If \code{dist}= "none", the mean is used as
     a fixed constant.}
  \item{btk}{list arguments for Bt/K, the relative depletive level in a specific reference year (\code{refyr}). 
   \code{dist} is the statistical distribution name from which to sample \code{btk}. 
   \code{low} and \code{up} are the lower and upper bounds of \code{btk} 
   in the selected distribution.  \code{mean} and \code{sd} are the mean and standard deviation 
    for selected distributions. The following are valid distributions: "none", "unif" - uniform, 
    "norm" - normal, "lnorm" - log-normal, "gamma" - gamma, and "beta" - beta distributions. 
    "unif" requires non-missing values for \code{low} and \code{up}. "norm", "lnorm", 
    "gamma" and "beta" require non-missing values for \code{low},\code{up}, \code{mean} and 
    \code{sd}. If "lnorm" is used, \code{mean} and \code{sd} must be on the natural log scale
    (keep \code{low} and \code{up} on the original scale). If \code{dist}= "none", the mean is used as
     a fixed constant. \code{refyr} is the selected terminal year  and can range from the first year 
     to the year after the last year of catch (t+1). 
  }
\item{fmsym}{list arguments for Fmsy/M. \code{dist} is the statistical distribution name from which 
   to sample \code{Fmsy/M}. \code{low} and \code{up} are the lower and upper bounds of \code{Fmsy/M} in
  the selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for selected 
  distributions. Valid distributions are the same as in \code{btk}. If \code{dist}= "none", the mean is used as
  a fixed constant. 
  }
\item{bmsyk}{list arguments for Bmsy/k. \code{dist} is the statistical distribution name from which 
   to sample \code{Bmsy/k}. \code{low} and \code{up} are the lower and upper bounds of \code{Bmsy/k} in
  the selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for selected 
  distributions. Valid distributions are the same as in \code{btk}. If \code{dist}= "none", the mean is used as
     a fixed constant. 
  }
\item{M}{list arguments for natural mortality. \code{dist} is the statistical distribution name from 
 which to sample \code{M}. \code{low} and \code{up} are the lower and upper bounds of \code{M} in the 
 selected distribution. \code{mean} and \code{sd} are the mean and standard deviation for selected 
 distributions. Valid distributions are the same as in \code{btk}. If \code{dist}= "none", the mean is used as
 a fixed constant. M is used to determine exploitation rate (Umsy) at MSY.
  }
\item{nsims}{number of Monte Carlos samples.} 
\item{catchout}{if catch is resampled, output the time series from every MC sample to a .csv file.
 0 = no (default), 1 = yes.} 
\item{grout}{numeric argument specifying whether graphs should be printed to console only (1) or to 
 both the console and TIF graph files (2).Use \code{setwd} before running function to direct .tif files
 to a specific directory. Each name of each file is automatically determined.} 
\item{graphs}{vector specifying which graphs should be produced. 1 = line plot of observed catch versus
  year, 2 = histogram of plausible (accepted) \code{k} values, 3 = histogram of plausible Bmsy values, 
 4 = histogram of plausible MSY values, 5 = histogram of plausible Fmsy values, 6 = histogram of Umsy values,
 7 = histogram of plausible Cmsy , 8 = histogram of Bmsy from plausible M, 9 = histogram of plausible Bt/k values,
 10 = histogram of plausible Fmsy/M values, 11 = histogram of plausible Bmsy/k values and 12 = histogram of
 plausible biomasses in \code{termyr}, 13 = line plots of accepted and rejected biomass trajectores
  with median and 2.5th and 97.5th percentiles (in red) and 14 =  stacked histograms of accepted and 
  rejected values for each input parameter and resulting estimates and if \code{grout}=2,
  .tif files are saved with "AR" suffix. Any combination of graphs can be 
  selected within c().  Default is all.
  }
\item{grargs}{list control arguments for plotting functions. \code{lwd} is the line width for graph = 1 and 13, 
 \code{nclasses} is the nclass argument for the histogram plots (graphs 2-12,14), \code{mains} and 
 \code{cex.main} are the titles and character expansion values for the graphs, \code{cex.axis} is the 
 character expansion value(s) for the x and y-axis tick labels and \code{cex.lab} is the character 
 expansion value(s) for the x and y-axis labels.  Single values of \code{nclasses},\code{mains}, 
 \code{cex.main},\code{cex.axis}, \code{cex.lab} are applied to all graphs.  To change arguments for 
 specific graphs, enclose arguments within c() in order of the number specified in \code{graphs}. 
}
\item{pstats}{list control arguments for plotting the median and 2.5% and 97.5% percentiles of parameters 
 and management quantities on respective graphs. \code{ol} = 0, do not overlay values on plots, 1 = 
 overlay values on plots. \code{mlty} and \code{mlwd} are the line type and line width of the median value;
  \code{llty} and \code{llwd} are the line type and line width of the 2.5% percentile; code{ulty} and 
\code{ulwd} are the line type and line width of the 97.5% percentile.
}
\item{grtif}{list arguments for the .TIF graph files. See \code{tiff} help file in R.} 
}
\details{
The method of Dick and MAcCall (2011) is used to produce estimates of MSY where only catch and 
information on resilience and current relative depletion is known. 

The  delay-difference model is used to propogate biomass: 
 
    B[t+1]<-B[t]+P[Bt-a]-C[t]

where \code{B[t]} is biomass in year t, \code{P[Bt-a]} is latent annual production based 
on parental biomass \code{agemat} years earlier and \code{C[t]} is the catch in year 
t. Biomass in the first year is assumed equal to \code{k}.

If Bmsy/k>=0.5, then P[t] is calculated as

P[t]<-g*MSY*(B[t-agemat]/k)-g*MSY*(B[t-agemat]/k)^n

where MSY is k*Bmsy/k*Umsy, n is solved iteratively using the equation, Bmsy/k=n^(1/(1-n)),
 and g is (n^(n/(n-1)))/(n-1). Fmsy is calculated as Fmsy=Fmsy/M*M and Umsy is calculated as 
(Fmsy/(Fmsy+M))*(1-exp(-Fmsy-M)).

If Bsmy/k < 0.5, Bjoin is calculated based on linear rules:
If Bmsy/k<0.3, Bjoin=0.5*Bmsy/k*k
If Bmsy/k>0.3 and Bmsy/k<0.5, Bjoin=(0.75*Bmsy/k-0.075)*k
 
If any B[t-a]<Bjoin, then the Schaefer model is used to calculated P:

P[Bt-agematt<Bjoin]<-B[t-agemat]*(P(Bjoin)/Bjoin+c(B[t-agemat]-Bjoin))

where c =(1-n)*g*MSY*Bjoin^(n-2)*K^(-n)
 
Biomass at MSY is calculated as: Bmsy=(Bmsy/k)*k

The overfishing limit (OFL) is Umsy*B[termyr].

\code{length(year)+1} biomass estimates are made for each run.

The rule for accepting a run is:
if(min(B)>0 && max(B)<=k && 

(objective function minimum<=tol^2) && abs(((max(B)-k)/k)*100)<=permax 
&& n<=maxn

If using the R Gui (not Rstudio), run 

      graphics.off()
  	windows(width=10, height=12,record=TRUE)
      .SavedPlots <- NULL

before running the dbsra function to recall plots.

}
\value{
  \item{Initial}{dataframe containing the descriptive statistics for each explored parameter.}
  \item{Parameters}{dataframe containing the mean, median, 2.5th and 97.5% percentiles, minumum and maximum
       of the plausible (accepted: likelihood(ll)=1) parameters. }  
  \item{Estimates}{dataframe containing the mean, median, 2.5th and 97.5% percentiles, minimum and maximum 
        of the management quantities (i.e., MSY, Bmsy, etc.) from the plausible parameters 
        (likelihood=1)}
  \item{Values}{dataframe containing the values of likelihood, k, Bt/k, Bmsy/k, M and associated management 
       quantities for all (likelihood=0 and likelihood=1) random draws.}
  \item{agemat}{agemat for use in function \code{dlproj}.}
 \item{end1yr}{value of the last year of catch data + 1 for use in function \code{dlproj}.}
 \item{type}{designates the output object as a \code{catchmsy} object for use in function \code{dlproj}.}

 The biomass estimates from each simulation are not stored in memory but are automatically 
     saved to a .csv file named "Biotraj-dbsra.csv". Yearly values for each simulation are stored across 
     columns.  The first column holds the likelihood values for each simulation (1= accepted, 0 = rejected).  
     The number of rows equals the number of simulations (\code{nsims}). This file is loaded to plot
      graph 13 and it must be present in the default or \code{setwd()} directory.  
    
  When \code{catchout}=1,   catch values randomly selected are saved to a .csv file named "Catchtraj-dbsra.csv". 
     Yearly values for each simulation are stored across columns.  The first column holds the likelihood 
     values (1= accepted, 0 = rejected).  The number of rows equals the number of simulations (\code{nsims}).   
  
 Use \code{setwd()} before running the function to change the directory where .csv files are stored.  
 
}
\references{

Dick, E. J. and A. D. MacCall. 2011. Depletion-based stock reduction analysis: a catch-based method for determining
sustainable yield for data-poor fish stocks. Fisheries Research 110: 331-341.
}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\note{The random distribution function was adapted from Nadarajah, S and S. Kotz. 2006. R 
 programs for computing truncated distributions. Journal of Statistical Software 16, code snippet 2.
}
\seealso{\code{\link{catchmsy}} \code{\link{dlproj}}}
\examples{
 \dontrun{
  data(cowcod)
  dbsra(year =cowcod$year, catch = cowcod$catch, catchCV = NULL, 
    catargs = list(dist="none",low=0,up=Inf,unit="MT"),
    agemat=11, k = list(low=100,up=15000,tol=0.01,permax=1000), 
    b1k = list(dist="none",low=0.01,up=0.99,mean=1,sd=0.1),
    btk = list(dist="beta",low=0.01,up=0.99,mean=0.1,sd=0.1,refyr=2009),
    fmsym = list(dist="lnorm",low=0.1,up=2,mean=-0.223,sd=0.2),
    bmsyk = list(dist="beta",low=0.05,up=0.95,mean=0.4,sd=0.05),
    M = list(dist="lnorm",low=0.001,up=1,mean=-2.90,sd=0.4),
    nsims = 10000)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
