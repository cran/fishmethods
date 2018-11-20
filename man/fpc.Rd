\name{fpc}
\alias{fpc}
\title{
Fishing Power Correction Factor from Experimental Fishing
}
\description{
Calculates fishing power correction ratios between two vessels or gears
}
\usage{
fpc(cpue1 = NULL, cpue2 = NULL, method = c(1,2,3,4),  deletezerosets = FALSE, 
kapp_zeros = "paired", boot_type = "paired", nboot = 1000, dint = c(1e-9,5),
 rint = c(1e-9, 20), decimals = 2, alpha = 0.05)
}
\arguments{
  \item{cpue1}{vector of CPUEs from vessel or gear considered the standard or baseline}
  \item{cpue2}{vector of CPUEs from other vessel or gear}
  \item{method}{method(s) to use to estimate fishing power correction. 1 = Ratio of Means, 2 = Randomized Block ANOVA,
  3 = Multiplicative Model, 4 = Kappenman 1992. Default = c(1,2,3,4)} 
  \item{deletezerosets}{if TRUE, paired observations with any CPUE=0 are eliminated prior to estimation. Default = FALSE.}
  \item{kapp_zeros}{for method = 4, how CPUE=0 is eliminated. "paired" eliminates the row of paired CPUE observations if CPUE = 0 is present
   for any observation within the pair,
  "ind" eliminates CPUE = 0 from  the individual CPUE vectors.}
  \item{boot_type}{the method for bootstrapping data. "paired" = resample paired CPUE observations, "unpaired" = resample individual CPUE vectors}
  \item{nboot}{the number of bootstrap replicates. Default = 1000.}
  \item{dint}{the lower and upper limits of the function interval searched by 
   function \code{uniroot} to solve Kappenman's \emph{d}.}
  \item{rint}{the lower and upper limits of the function interval searched by 
   function \code{optimize} to solve Kappenman's \emph{r}.}
  \item{decimals}{the number of decimal places for output of estimates.}
  \item{alpha}{the alpha level used to calculate confidence intervals.}
}
\details{
 
The four methods for estimating fishing power correction factors given in Wilderbuer et al. (1998) are encoded.

If paired CPUE observations are both zero, the row is automatically eliminated. If \code{deletezerosets} = TRUE, the paired
CPUE observations with any CPUE = 0 will be eliminated. 

Zeroes are allowed in methods 1, 2 and 3. 

For the Kappenman method (method=4), only non-zero CPUEs are allowed. Use \code{kapp_zeros} to select the elimination
 method. An unequal number of observations between vessels is allowed in this method and can result using 
\code{kapp_zeros} = "ind".  FPC is derived by using the methodology where r that minimizes the sum of
squares under the first conjecture relative to the second is estimated (Kappenman 1992: 2989; von Szalay and Brown 2001). 


Standard errors and confidence intervals of FPC estimates are derived for most methods by using an approximation formula (where applicable),
jackknifing and/or bootstrapping. Specify the type of bootstrapping through \code{boot_type}. For methods 1-3,
jackknife estimates are provided only when \code{boot_type}="paired". If method = 4, jackknife estimates are provided only
when \code{boot_type}="paired" and \code{kapp_zeros}="paired". 

Confidence intervals are provided for the approximation formulae specified in Wilderbuer et al (1998), the jackknife estimates and bootstrap
estimates. Confidence intervals for the jackknife method are calculated using the standard formula
(estimate+/-z[alpha/2]*jackknife standard error).  Bootstrap confidence intervals are derived using the percentile method 
(Haddon 2001).

}
\value{
A dataframe containing method name, sample size for cpue1 (n1) and cpue2 (n2) ,mean cpue1,
mean cpue2, fishing power correction (FPC), standard error from approximation formulae (U_SE),
standard error from jackknifing (Jack_SE), standard error from bootstrapping (Boot_SE), lower and upper confidence intervals
from approximation formulae (\code{U_X\%_LCI} and \code{U_X\%_UCI}),lower and upper confidence intervals
from jackknifing (\code{Jack_X\%_LCI} and \code{Jack_X\%_UCI}) and lower and upper confidence intervals
from bootstrapping (\code{Boot_X\%_LCI} and \code{Boot_X\%_UCI}).
}
\references{

Haddon, M. 2001. Modelling and Quantitative Methods in Fisheries. Chapman & Hall/CRC Press. Boca Raton, Florida.


Kappenman, R. F. 1992. Robust estimation of the ratio of scale parameters for positive
random variables. Communications in Statistics, Theory and Methods 21: 2983-2996.


von Szalay, P. G. and E. Brown. 2001. Trawl comparisons of fishing power differences and
their applicability to National Marine Fisheries Service and Alask Department of Fish
and Game trawl survey gear.  Alaska Fishery Research Bulletin 8(2):85-95.


Wilderbuer, T. K., R. F. Kappenman and D. R. Gunderson. 1998. Analysis of fishing power correction factor 
estimates from a trawl comparison experiment.  North American Journal of Fisheries Management 18:11-18.

}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\examples{
\dontrun{
 #FPC for flathead sole from von Szalay and Brown 2001
   data(sole)
   fpc(cpue1=sole$nmfs,cpue2=sole$adfg,boot_type="unpaired",kapp_zeros="ind",method=c(4),
            alpha=0.05)}
}
\keyword{misc}
