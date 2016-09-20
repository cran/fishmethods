\name{fm_telemetry}
\alias{fm_telemetry}
\title{Estimation of Fishing and Natural Mortality from Telemetry Data
}
\description{
The method of Hightower et al. (2001) is implemented to estimate fishing mortality, natural mortality and probability of detection from telemetry data.
}
\usage{
fm_telemetry(filetype = c(1), caphistory = NULL, Fdesign = NULL, Mdesign = NULL, 
Pdesign = NULL, whichlivecells =  NULL, 
whichdeadcells = NULL, constant = 1e-14, initial = NULL, 
invtol = 1e-44, control = list(reltol=1e-8,maxit=1000000))
}
\arguments{
\item{filetype}{type of file to read. 1 = R character vector with individual capture histories (1 history per row), 
or 2 = an external text file with individual capture histories. If \code{filetype}=2, then the capture histories 
in the file should not be enclosed in quotes and there 
should not be a column name. 
}
\item{caphistory}{File or R object with capture histories. If filetype=2, location and filename of text file enclosed in quotes (e.g., \dQuote{C:/temp/data.txt}).
}
\item{Fdesign}{vector of characters specifying the occasion parameter structure for fishing mortality (F). See details.
}
\item{Mdesign}{vector of characters specifying the occasion parameter structure for natural mortality (M). See details.
}
\item{Pdesign}{vector of characters specifying the occasion parameter structure for the probability of detection (P). See details.}
\item{whichlivecells}{list containing the structure of occasion live cells to use in each release during the estimation process. Multiple ranges may be specified. 
For each range, specify the first release, last release, and number of observed occasions (cells) enclosed within \code{c()}. For example, to use the first 4 cells of releases 1-5, specify \code{c(1,5,4)}. \code{whichlivecells} is a list object of all ranges (e.g., whichlivecells =list(c(1,5,4),c(6,26,6))). Specify \code{whichlivecells=NULL} to use all cells. 
The Hightower et al. (2001) specification is 
\code{whichlivecells}=list(c(1,5,4),c(6,6,5),
c(7,26,4)).
} 
\item{whichdeadcells}{list containing the structure of occasion dead cells to used in each release during the estimation process. Same as \code{whichlivecells}. The Hightower et al. (2001) specification is 
\code{whichdeadcells}=list(c(1,5,4),c(6,6,6),
c(7,26,4))
}
\item{constant}{A small number to use in the multinomial log-likelihood (Obs * 
log(max(constant,
Expected Prob))) to avoid errors if any probability is 0. If the number is too large, it may affect the minimization of the likelihood. Default is 1e-14.
} 
\item{initial}{vector of starting values for fishing and natural mortality, and the probability of detection. First position is the starting value for all Fs, the second position is the starting value for all Ms, and the third position is the starting value for all Ps (e.g., c(0.1,0.2,0.8)).
}
\item{invtol}{the tolerance for detecting linear dependencies in the columns of a in \code{solve}(the function used to invert the hessian matrix). Adjust this value if errors about tolerance limits arise.
}
\item{control}{A list of control parameters for \code{optim}.     See function \code{optim} for details.}
}
\details{
The telemetry method of Hightower et al. (2001) is implemented. Individual capture histories are used in the function.  The function uses complete capture histories (Burnham et al., 1987) and it is the presence of specific codes in the individual capture histories that split the capture histories into live and dead arrays. F and M estimates are needed for occasions 1 to the total number of occasions minus 1 and P estimates are needed for occasions 2 to the total number of occasions.

Capture histories are coded following Burnham et al. (1987)(i.e., 0 = not relocated, and 1 = relocated) with the following exceptions:

All live relocations are coded with 1. If a fish is relocated and is dead, then \code{D} is used. For example,

\code{101011}  - fish released on occasion 1 is relocated alive on occasions 3,5 and 6

\code{10111D}  - fish released on occasion 1 is relocated alive on occasions 3,4,and 5 but is relocated dead on occasion 6.

New releases are allowed to occur on multiple occasions.  The capture history of newly-released individuals should be coded with a zero (0) for the occasions before their release.

\code{100110} - fish released on occasion 1 is relocated live on occasion 4 and 5 

\code{101000} - fish released on occasion 1 is relocated live on occasion 3

\code{010111} - fish released on occasion 2 is relocated live on occasion 4, 5 and 6

\code{011000} - fish released on occasion 2 is relocated live on occasion 3

\code{001101} - fish released on occasion 3 is relocated live on occasion 4 and 6

\code{00100D} - fish released on occasion 3 is relocated dead on occasion 6.

To censor fish from the analyses, specify \code{E} after the last live encounter. For example,

\code{10111E000}  - fish released on occasion 1 is relocated alive on occasions 3,4,and 5 but is believed to have emigrated from the area by occasion 6. The capture history before the \code{E} will be used, but the fish is not included in the virtual release in occasion 6. 

All life histories are summarized to reduced m-arrays (Burnham et al. (1987: page 47, Table 1.15). 

The function \code{optim} is used to find F, M and P parameters that minimize the negative log-likelihood. Only cells specified in \code{whichlivecells} and \code{whichdeadcells} are used in parameter estimation.

The logit transformation is used in the estimation process to constrain values between 0 and 1. Logit-scale estimated parameters are used to calculate Sf=1/(1+exp(-B)), Sm=1/(1+exp(-C)) and P=1/(1+exp-(D)).  F and M are obtained by -log(Sf) and -log(Sm).

The standard error of Sfs, Sm, P, F and M are obtained by the delta method: 

SE(Sf)=sqrt((var(B)*exp(2*B))/(1+exp(B))^4), 

SE(Sm)=sqrt((var(C)*exp(2*C))/(1+exp(C))^4),

SE(P)=sqrt((var(D)*exp(2*D))/(1+exp(D))^4),

SE(F)=sqrt(SE(Sf)^2/Sf^2),

SE(M)=sqrt(SE(Sm)^2/Sm^2).

All summary statistics follow Burnham and Anderson (2002).  Model degrees of freedom are calculated as nlive+ndead+nnever-nreleases-1-npar where nlive is the number of \code{whichlivecells} cells, ndead is the number of \code{whichdeadcells} cells, nnever is the number of never-seen cells, nreleases is the number of releases and npar is the number of estimated parameters. Total chi-square is calculated by summing the cell chi-square values.

The program allows the configuration of different model structures (biological realistic models) for the estimation of fishing and natural mortalities, and detection probabilities. These structures are specified in \code{Fdesign}, \code{Mdesign} and \code{Pdesign}. Consider the following examples:

\emph{Example 1}

Tags are relocated over seven occasions.  One model structure might be constant fishing mortality estimates over occasions 1-3 and 4-6, one constant estimate of natural mortality for the entire sampling period, and one estimate of probability of detection for each occasion. To specify this model structure:
\code{Fdesign} is c(\dQuote{1},\dQuote{4}), \code{Mdesign} is c(\dQuote{1}) and the \code{Pdesign} is c(\dQuote{2:2}).

Note: The structures of \code{Fdesign} and \code{Mdesign} must always start with the first occasion, whereas the structure for \code{Pdesign} must always start with the second occasion. 

Use the multiplication sign to specify occasions whose estimates of F, M or P will be taken from values of other occasions.

\emph{Example 2}

Tags are relocated over six occasions.  One model structure might be separate fishing mortality estimates for occasions 1-3 but assign the same parameter estimates to occasions 4-6, one constant estimate of natural mortality for occasions 1-5 and 6, and one constant probability of detection over all occasions.   The \code{Fdesign} is c(\dQuote{1:3*4:6}), the \code{Mdesign} is c(\dQuote{1},\dQuote{6}) and the \code{Pdesign} is c(\dQuote{2}).

\emph{Example 3}

Specification of model 18 listed in Table 1 of Hightower et al. (2001) is shown.  Each occasion represented a quarter of the year. The quarter-year design for F, M and P specifies that quarterly estimates are made in each year. \code{Fdesign} is c(\dQuote{1},\dQuote{4},\dQuote{7},\dQuote{11},\dQuote{14},\dQuote{17}, 
\dQuote{20},\dQuote{24}). \code{Mdesign} is c(\dQuote{1},\dQuote{4},\dQuote{7},\dQuote{11},\dQuote{14},\dQuote{17},\dQuote{20},\dQuote{24}) and the \code{Pdesign} is c(\dQuote{2},\dQuote{4},\dQuote{7},\dQuote{11},\dQuote{14},\dQuote{17}, 
\dQuote{20},\dQuote{24}). 

If the number of occasions to be assigned parameters from other occasions are less than the number of original parameters (e.g., c(\dQuote{11:13*24:25}), then only the beginning sequence of original parameters equal to the number of occasions are used. For instance, in c(\dQuote{11:13*24:25}), only parameters 11 and 12 would be assigned to occasions 24 and 25.

If the number of occasions to be assigned parameters from other occasions are greater than the number of original parameters (e.g., c(\dQuote{11:12*24:26})), then the last original parameter is re-cycled. In the example c(\dQuote{11:12*24:26}), the parameter for occasion 12 is assigned to occasions 25 \emph{and} 26.

To assist with the parameter structures, function \code{fm_checkdesign} may be used to check the desired design before use in this function. 

If values of standard error are NA in the output, the hessian matrix used to claculate the variance-covariance matrix could not be inverted. If this occurs, try adjusting the \code{reltol} argument (for more information, see function \code{optim}).

In this function, the never-seen expected number is calculated by summing the live and dead probabilities, subtracting the number from 1, and then multiplying it by the number of releases. No rounding occurs in this function.
 
The multinomial likelihood includes the binomial coefficient.

Model averaging of model can be accomplished using the function \code{fm_model_avg}.

Note: In Hightower et al.'s original analysis, the cell probability code in SURVIV for the dead relocation in release occasion 6 had an error. The corrected analysis changed the estimates for occasions 11-13 compared to the original published values.
}
\value{List containing summary statistics for the model fit, model convergence status, parameter estimates estimates of fishing mortality, natural mortality, and probabilties of detection and standard errors by occasion, the parameter structure (Fdeisgn, Mdesign and Pdesign), the m-arrays, the expected (predicted) number of live and dead relocations,  cell chi-square and Pearson values for live and dead relocations,  matrices with the probability of being relocated alive and dead by occasion, the whichlivecells and whichdeadcells structures, and configuration label (type) used in the \code{fm_model_avg} function.
}
\references{ 
Burnham, K. P. and D. R. Anderson. 2002. Model selection and multimodel inference : A Practical Information-Theorectic
 Approach, 2nd edition. Spriner-Verlag, New York, NY. 488 p.

Burnham, K. P. D. R. Anderson, G. C. White, C. Brownie, and K. H. Pollock. 1987. Design and analysis methods for fish survival experiments based on release-recapture.  American FIsheries Society Monograph 5, Bethesda, Maryland.

Hightower, J. E., J. R. Jackson, and K. H. Pollock. 2001. Use of telemetry methods to estimate natural and fishing mortality of striped bass in Lake Gaston, North Carolina. Transactions of the American Fisheries Society 130: 557-567.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}
}
\seealso{\code{fm_model_avg},\code{fm_checkdesign}}
\examples{
\dontrun{
# Set up for Full model of Hightower et al.(2001)
data(Hightower)
fm_telemetry(filetype=1,caphistory=Hightower$caphistory, Fdesign=c("1:26"),
 Mdesign=c("1:26"),Pdesign = c("2:25"),
whichlivecells=list(c(1,5,4), c(6,6,5),
 c(7,26,4)), 
whichdeadcells=list(c(1,5,4), c(6,6,6),
 c(7,26,4)),
initial=c(0.05,0.02,0.8),
control=list(reltol=1e-5,maxit=1000000))

#Set up for best model F(Qtr,yr), M constant, Pocc
fm_telemetry(filetype=1,caphistory=Hightower$caphistory, Fdesign=c("1", "4", "7", "11",
 "14", "17", "20", "24"), 
Mdesign=c("1"), Pdesign = c("2:27"),
whichlivecells=list(c(1,5,4), c(6,6,5),
 c(7,26,4)),
whichdeadcells=list(c(1,5,4), c(6,6,6),
 c(7,26,4)), 
initial=c(0.05,0.02,0.8),
 control=list(reltol=1e-8,maxit=1000000))
}
}            
\keyword{misc}
