\name{fm_checkdesign}
\alias{fm_checkdesign}
\title{Check parameter structure of Hightower et al. (2001) models}
\description{
Check design of parameter structure before use in function \code{fm_telemetry}.
}
\usage{
fm_checkdesign(occasions = NULL, design = NULL, type = "F" )
}
\arguments{
  \item{occasions}{total number of occasions that will be modeled in data}
  \item{design}{vector of characters specifying the occasion parameter structure (see details).}
  \item{type}{character type of parameter to which design will be applied: F = fishing mortality, M = natural mortality, and P = probability of detection. Default = F.}
}
\details{The program allows the configuration of different parameter structure for the estimation of fishing and natural mortalities, and detection probabilities. These structures are specified in \code{design}.  Consider the following examples:

\emph{Example 1}

Tags are relocated over seven occasions.  One model structure might be constant fishing mortality estimates over occasions 1-3 and 4-6. To specify this model structure:
\code{design} is c(\dQuote{1},\dQuote{4}).

Note: The structures of \code{design} must always contain the first occasion for fishing mortality and natural mortality, whereas the structure for the probability of detection must not contain the first occasion. 

\emph{Example 2}

Tags are relocated over six occasions.  One model structure might be separate fishing mortality estimates for occasion 1-3 and the same parameter estimates for occasions 4-6. The \code{design} is c(\dQuote{1:3*4:6}).

Note: The structures of \code{Fdesign} and \code{Mdesign} must always start with the first occasion, whereas the structure for \code{Pdesign} must always start with the second occasion. 

Use the multiplication sign to specify occasions whose estimates of F, M or P will be taken from values of other occasions.

\emph{Example 3}

Specification of model 3 listed in Table 1 of Hightower et al. (2001) is shown.  Each occasion represented a quarter of the year. The quarter design for F specifies that quarterly estimates are the same in both years. \code{design} is c(\dQuote{1*14},\dQuote{4*17},\dQuote{7*20},\dQuote{11*24}).

\emph{Example 4}

In Hightower et al. (2001), the quarter and year design specifies that estimates are made for each quarter but are different for each year. \code{design} is 

c(\dQuote{1}, \dQuote{4}, \dQuote{7}, \dQuote{11}, 
\dQuote{14}, \dQuote{17}, \dQuote{20}, \dQuote{24}).

If the number of occasions to be assigned parameters from other occasions are less than the number of original parameters (e.g., c(\dQuote{11:13*24:25}), then only the beginning sequence of original parameters equal to the number of occasions are used. For instance, in c(\dQuote{11:13*24:25}), only parameters 11 and 12 would be assigned to occasions 24 and 25.

If the number of occasions to be assigned parameters from other occasions are greater than the number of original parameters (e.g., c(\dQuote{11:12*24:26})), then the last original parameter is re-cycled. In the example c(\dQuote{11:12*24:26}), the parameter for occasion 12 is assigned to occasions 25 \emph{and} 26.
}
\value{dataframe containing the parameter order by occasion. 
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{fm_telemetry}}}
\examples{
fm_checkdesign(occasions=27, design=c("1*14","4*17","7*20","11*24"),type="F")
}            
\keyword{misc}

