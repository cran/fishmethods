\name{LengthIncr}
\alias{LengthIncr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Re-parameterized Von Bertalanffy growth function (GROTAG) from Francis (1988)
developed in the context of analysis of tagging data
}
\description{
This function calculates the expected growth increment of a fish of
length L1 over a time period of delta.T, assuming the knowledge of the
mean annual growth rates (g.alpha and g.beta) of fish of arbitrary length alpha and beta.
}
\usage{
LengthIncr(gr.alpha, gr.beta, L1, delta.T, alpha = 35, beta = 55)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr.alpha}{
A numeric value giving the mean annual growth rate of fish of arbitrary
length alpha 
}
  \item{gr.beta}{
A numeric value giving the mean annual growth rate of fish of arbitrary
length beta 
}
  \item{L1}{
A vector giving the length at release of a tagged fish
}
  \item{delta.T}{
A vector giving the time at liberty of a fish after tagging
}
  \item{alpha}{
A numeric value giving an arbitrary length alpha
}
  \item{beta}{
A numeric value giving an arbitrary length beta (NB: beta > alpha)
}
}
\details{
 Use this function to fit the models of Francis (1988).  Examples of how this is accomplished are given below in the EXAMPLES section.
}
\references{
Francis, R.I.C.C., 1988. Maximum likelihood estimation of growth and growth variability from tagging data. New Zealand Journal of Marine and Freshwater Research, 22, p.42--51. 
}
\author{
Marco.Kienzle@gmail.com
}
\note{
It is assumed that the arbitrary length alpha > beta
}
\examples{
data(bonito)

##### Example of implementation of Francis (1988) model 1

# Plot the data
plot(bonito$L1, bonito$delta.L,xlab = "Length at tagging", ylab = "Length increment", las=1, main = "Bonito data from Campbell et al. (1975)")

# function
f1 <-  function(p){
sum(dnorm(x=bonito$delta.L, mean = LengthIncr(p[1], p[2], bonito$L1, delta.T = bonito$T2 - bonito$T1, alpha = 35, beta = 55), sd = p[3], log = TRUE))
}

# Get initial values for g35 and g55
# by calculating average growth from data with T2-T1 superior to 0.1 year
g35.init =  with(subset(bonito, L1 > 30 & L1 < 40 & (T2-T1) > 0.1),  mean(delta.L/ (T2 - T1)))
g55.init =  with(subset(bonito, L1 > 50 & L1 < 60 & (T2-T1) > 0.1),  mean(delta.L / (T2 - T1)))

# Fit the model - maximum likelihood parameters estimation
mod1 = optim(par = c(g35 = g35.init, g55 = g55.init, s = 0.5), fn = f1,
lower = c(max(g35.init - 5, 0), max(g55.init - 5,0), 0), upper =
c(g35.init + 5, g55.init + 5, 5), method = "L-BFGS-B", hessian = TRUE,
control = list(fnscale = -1))



##### Example of implementation of Francis (1988) model 2
##### model with g.alpha, g.beta, standard deviation of error measurement s and growth variability as a function of mean growth

# function - Assume standard deviation of delta.L is proportional to length increment but replace negative values the value for measurement error
f2 <- function(p){
std.dev = LengthIncr(p[1], p[2], bonito$L1, delta.T = bonito$T2 - bonito$T1, alpha = 35, beta = 55)
std.dev = replace(std.dev, which(std.dev < 0), p[4])
sum(dnorm(x = bonito$delta.L, mean = LengthIncr(p[1], p[2], bonito$L1, delta.T = bonito$T2 - bonito$T1, alpha = 35, beta = 55 ), sd = sqrt( (p[3] * std.dev) ^ 2 + p[4] ^ 2), log = TRUE))
}

# Get first hint on g35 and g55 by estimating the average growth using the data that have T2-T1 superior to 0.1 year
g35.init =  with(subset(bonito, L1 > 30 & L1 < 40 & (T2-T1) > 0.1),  mean(delta.L/ (T2 - T1)))
g55.init =  with(subset(bonito, L1 > 50 & L1 < 60 & (T2-T1) > 0.1),  mean(delta.L / (T2 - T1)))

# Fit the model - maximum likelihood parameters estimation
mod2 = optim(c(g35 = g35.init, g55 = g55.init, factor = 0.5, s = 0.5), f2, lower = c(max(g35.init - 5, 0), max(g55.init - 5,0), 0, 0), upper = c(g35.init + 5, g55.init + 5, 5, 5), method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1))



##### Example of implementation of Francis (1988) model 3
##### model with g.alpha, g.beta, mean and standard deviation of error measurement s and growth variability as a function of mean growth

# function - Assume standard deviation of delta.L is proportional to length increment but replace negative values the value for measurement error
f3 = function(p){
std.dev = LengthIncr(p[1], p[2], bonito$L1, delta.T = bonito$T2 - bonito$T1, alpha = 35, beta = 55)
std.dev = replace(std.dev, which(std.dev < 0),p[5])
sum(dnorm(x = bonito$delta.L, mean = LengthIncr(p[1], p[2], bonito$L1, delta.T = bonito$T2 - bonito$T1, alpha = 35, beta = 55) + p[4] , sd = sqrt( (p[3] *std.dev) ^ 2 + p[5] ^ 2), log = TRUE))
}

# Get first hint on g35 and g55 by estimating the average growth using the data that have T2-T1 superior to 0.1 year
g35.init =  with(subset(bonito, L1 > 30 & L1 < 40 & (T2-T1) > 0.1),  mean(delta.L/ (T2 - T1)))
g55.init =  with(subset(bonito, L1 > 50 & L1 < 60 & (T2-T1) > 0.1),  mean(delta.L / (T2 - T1)))

# Fit the model - maximum likelihood parameters estimation
mod3 = optim(c(g35 = g35.init, g55 = g55.init, factor = 1, m = -1, s = 0.5), f3, lower = c(max(g35.init - 5, 0), max(g55.init - 5,0), 0, -2,0), upper = c(g35.init + 5, g55.init + 5, 5, 2,5), method = "L-BFGS-B", hessian = TRUE, control = list(fnscale = -1))



##### Example of implementation of Francis (1988) model 5
##### Model with 4 parameters: g.alpha, g.beta, measurement error and outliers

# function 
f5 = function(p){
# Calculate range of growth increment
R = range(bonito$delta.L)[2] - range(bonito$delta.L)[1] 
# Likelihood function
sum(log((1-p[3]) * dnorm(x = bonito$delta.L, mean = LengthIncr(p[1], p[2], bonito$L1, delta.T = bonito$T2 - bonito$T1, alpha = 35, beta = 55), sd = p[4]) + p[3] / R)) 
}

# Get first hint on g35 and g55 by estimating the average growth using the data that have T2-T1 superior to 0.1 year
g35.init =  with(subset(bonito, L1 > 30 & L1 < 40 & (T2-T1) > 0.1),  mean(delta.L/ (T2 - T1)))
g55.init =  with(subset(bonito, L1 > 50 & L1 < 60 & (T2-T1) > 0.1),  mean(delta.L / (T2 - T1)))

# Fit the model - maximum likelihood parameters estimation
mod5 = optim(c(g35 = g35.init, g55 = g55.init, p = 0, s = 2), f5, lower
= c(max(g35.init - 5, 0), max(g55.init - 5,0),0,0), upper = c(g35.init +
5, g55.init + 5,1,5), method = "L-BFGS-B", hessian = TRUE, control =
list(fnscale = -1))



##### Generate a table formatted in the same way as Francis (1988) Tab. 3, p.47

# Header
result = as.data.frame(matrix(NA, ncol=1, nrow = 9))
dimnames(result)[[1]] = c("Log likelihood","Mean growth rates g35","Mean growth rates g55", "Seasonal variation u", "Seasonal variation v", "Growth variability, nu", "Measurement error, s", "Measurement error, m", "Outliers, p")

# mod1
result = cbind(result[,-1], model1 = c(round(mod1$value,1), round(mod1$par[1],1), round(mod1$par[2],1),NA,NA,NA,round(mod1$par[3],2),NA,NA))

# mod2
result = cbind(result, model2 = c(round(mod2$value,1), round(mod2$par[1],1), round(mod2$par[2],1),NA,NA,round(mod2$par[3],2),round(mod2$par[4],2),NA,NA))

# mod3
result = cbind(result, model3 = c(round(mod3$value,1), round(mod3$par[1],1), round(mod3$par[2],1),NA,NA,round(mod3$par[3],2),round(mod3$par[5],2),round(mod3$par[4],2),NA))

# mod5
result = cbind(result, model5 = c(round(mod5$value,1), round(mod5$par[1],1), round(mod5$par[2],1),NA,NA,NA,round(mod5$par[4],2),NA, round(mod5$par[3],3)))
print(result)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

