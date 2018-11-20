\name{fm_model_avg}
\alias{fm_model_avg}
\title{Model Averaging for the Telemetry Method of Hightower et al. (2001)}
\description{
 Calculates model averaged estimates of instantaneous fishing, natural and probability of detection for telemetry models of Hightower et al. (2001).}
\usage{
fm_model_avg(..., global = NULL, chat = 1)
}
\arguments{
  \item{\dots}{model object names separated by commas}
  \item{global}{specify global model name in quotes.  If the global model is the first model included 
     in the list of candidate models, this argument can be ignored.}
  \item{chat}{chat for the global model.}

}
\details{Model estimates are generated from function \code{fm_telemetry}.
Averaging of model estimates follows the procedures in Burnham and Anderson (2002). Variances of parameters are adjusted for overdispersion using the c-hat estimate from the global model 
: \code{sqrt(var*c-hat)}.  If c-hat of the global model is <1, then c-hat is set to 1. The c-hat is used to calculate the quasi-likelihood AIC and AICc 
metrics for each model (see page 69 in Burnham and Anderson(2002)). QAICc differences among models are calculated by
subtracting the QAICc of each model from the model with the smallest QAICc value. These differences are used to calculate 
the Akaike weights for each model following the formula on page 75 of Burnham and Anderson (2002). The Akaike weights are
used to calculate the weighted average and standard error of parameter estimates by summing the product of the model-specific Akaike weight and parameter estimate 
across all models.  An unconditional standard error is also calculated by 
\code{sqrt(sum(QAICc wgt of model i 
* (var of est of model i 
+ (est of model i - avg of all est)^2)))}.
}
\value{List containing model summary statistics, model-averaged estimates of fishing, natural and probability of detections and their weighted and uncondtional standard errors . 
}
\references{ 
Burnham, K. P. and D. R. Anderson. 2002. Model selection and multimodel inference : A Practical Information-Theorectic Approach, 2nd edition. Spriner-Verlag, New York, NY. 488 p.
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\seealso{\code{\link{fm_telemetry}}}
\examples{
## This is a typical specification, not a working example
\dontrun{
fm_model_avg(model1,model2,model3,model4,model5,model6,model7,global="model7")}
}            
\keyword{misc}

