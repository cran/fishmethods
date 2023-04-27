\name{compare.lrt.plus}
\alias{compare.lrt.plus}
\title{Comparison of growthlrt.plus model objects}
\description{
 Compute likelihood ratio tests for two or more growthlrt.plus model objects via Kimura (1990)}
\usage{
compare.lrt.plus(...)
}
\arguments{
  \item{\dots}{growthlrt.plus object names separated by commas}
}
\details{
Likelihood ratio and F tests are computed for models compared against one another in the order specified.
}
\value{List containing model test statistics 
}
\references{ 
Kimura, D. K. 1990. Testing nonlinear reression parameters under heteroscedastic, normally distributed errors.
Biometrics 46: 697-708.

}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{growthlrt.plus}}}
\examples{
## This is a typical specification, not a working example
\dontrun{
compare.lrt.plus(model1,model2)}
}            
\keyword{misc}

