\name{growthlrt.plus}
\alias{growthlrt.plus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Likelihood Methods for Comparing Multiple Growth Curves
}
\description{
Additional likelihood methods for comparison of two or more curves under heteroscedastic,
 normally-distributed errors and 
differing within-group variances based on Kimura (1990).
}
\usage{
growthlrt.plus(model, data, params = NULL, start = NULL, within_grp_var = ~1,
      cfh = NULL, nlminb.control = list(iter.max = 10000, rel.tol = 1e-10),
      optim.control=list(maxit = 10000, reltol = 1e-10))
}
\arguments{
  \item{model}{a two-sided formula object describing the model, with the response on the 
	left of a ~ operator and a nonlinear expression involving parameters on the right.}
  \item{data}{A data frame containing the variables named in \emph{model}. Rows should represent individual observations.}
  \item{params}{a two-sided linear formula of the form \code{p1=~1} or \code{p1=~group} for each parameter estimated in model. 
     The \code{p1} represents a parameter included on the right hand side of model. A \code{1} on the right hand 
	side of the formula indicates a single parameter is estimated, whereas a variable name of a group variable 
      will estimate as many parameters as there are levels in the group variable.}
  \item{start}{a required named list with the initial values for the parameters in model. If multiple
     estimates for a given parameter are desired, starting values should be enclosed in \code{c()}.}
  \item{within_grp_var}{a one-sided formula of the form \code{within_grp_var= ~1} or    
	\code{within_grp_var= ~group}. A \code{1} on the right hand side of the formula indicates a 
        single within-group variance is estimated for all groups, whereas a variable name 
        (same one used in \code{params}) will estimate different sigmas for each level under group.
	}
  \item{cfh}{NULL or a named list with arguments needed to correct for heterogeneity of variances. 
	If the latter, the required arguments are \code{form}, \code{value},and \code{fixed}. See details for more information.}
  \item{nlminb.control}{Additional controls passed to the optimization function \emph{nlminb}.}
  \item{optim.control}{Additional controls passed to the optimization function \emph{optim}.}
}
\details{
The likelihood methods of Kimura (1990) are used to fit any nonlinear equation, correct for 
heterogeneity of variances, and estimate common or separate within-group variances depending on
user-specifications. A main assumption is errors are normally-distributed.  
The results of the model fits can then be used in function \emph{compare.lrt.plus} 
to determine if parameterizations differ significantly from each other by using a likelihood ratio 
and an F test.

Steps of the modeling process are as follows:

1) Specify the nonlinear model equation in the same 
formula format as would be done in function \emph{nls}. For example, the von Bertalanffy growth equation 
is written as:

\code{sl~Linf*(1-exp(-K*(age-t0)))}

where \emph{sl} is the variable name for length data, \emph{age} is the variable name for age data,
and \emph{Linf}, \emph{K} and \emph{t0} are parameters to be estimated. 

2) Specify the parameter formulae under \code{params}.  These formulae 
are used to indicate that additional parameters based on some group variable should be estimated. 
For example,

\code{params=list(Linf~1,K~1,t0~1)}

specifies single parameters are estimated for Linf, K and t0.  

\code{params=list(Linf~sex,K~1,t0~1)}

specifies that separate Linfs are to be estimated for each sex and only single estimates
for K and t0.
 
\code{params=list(Linf~sex,K~sex,t0~sex)}

specifies that separate Linfs, Ks and t0s are to be estimated for each sex.  Different group variables 
for each parameter are not allowed.

3) Specify start values for all parameters.  For example, if separate Linfs, Ks and t0s for a group
variable like sex (only two-levels: M and F), then 6 starting values must be given.  When parameters are
based on a group variable, then the first estimate of a parameter will be the reference value (labeled as Intercept)
and the remaining parameters will be estimated as a deviation from that reference value. Reference values 
are determined by alphanumeric order of levels within the group variable.

\code{start=list(Linf=c(300,10),K=c(0.3,0.05),t0=c(0,-0.5))}

is an example of the starting values for the 6 parameter model mentioned above. Warning messages are generated
if the number of start values is less than or greater than the number of parameters being estimated. Internally,
code will add (1/10th of first value) or truncate (last number(s) in list) start values in these cases. 
However, the user should specify the appropriate number of values to ensure successful optimization. 

4) Specify the within-group variance structure.  If the within-group variance is assumed
the same among groups, then 

\code{within_grp_var=~1} 

which is the default specification. If within-group variances are suspected to differ among groups (e.g., sex), then

\code{within_grp_var=~sex}

Separate variances will be estimated for each level of the group variable.  Whether or not better model fits can be 
obtained by estimating separate group variances can be tested using the model comparison methods (see below).
When estimating thetas (correcting for heterogeneity), explore different starting values for the main parameters to 
ensure global convergence.



5) Specify the correction for heterogenity (\code{cfh}) argument(s) if needed.  Initial curve fittings should be performed 
and plots of residuals versus fitted values examined to determine if there is a change in residual variation with 
increasing fitted values. If so, this indicates the presense of heterogeneity in variance which must be corrected to obtain 
unbiased parameters estimates, standard errors, residual sum-of-squares, etc. Kimura (1990) uses the power function (same as the
 \emph{varPower} function in Pinheiro and Bates (2004)) and additional parameters known as \emph{theta} are estimated.
If \code{cfh} is NULL, then homogeneity of variance is assumed. If heterogeneity of variance needs to be accounted for, specify 
\code{cfh} as

\code{cfh=list(form=~1,value=0,fixed=NULL)}

\code{form} is a formula and is 1 if a single theta is assumed equal among groups. If individual thetas are desired for groups (heterogenity is different
for each group), then a group variable is used (e.g.,\emph{form}=~sex). 

\code{value} is the initial starting value(s) for theta(s). If more than 1 theta will
be estimated, provide the same number of starting values  within \emph{c()} as thetas.

\code{fixed} is used to indicate whether the thetas will be estimated (default \emph{NULL}) or assumed fixed to numeric values specified by the user.

\code{cfh=list(form=~sex,value=0,fixed=c(0.5,0.9))}

indicates that thetas for each sex (two-levels: M and F) will not be estimated, but fixed to values of 0.5 and 0.9
 

6) Run the model function.  Parameter estimation is performed intially by using the optimization function \emph{nlminb}. 
The estimated parameters are then used as starting values and optimization is performed again by 
using function \emph{optim} to obtain the final parameter estimates and the Hessian matrix from 
which standard errors are derived. Unlike estimation of thetas conducted in function \emph{gnls} in package \emph{nlme}, thetas herein
are estimated as parameters, standard errors are generated, and t-tests for significance are conducted. These extra parameters are counted in the
determination of residual and model degrees of freedom.   

To convert a non-reference level estimate to the same scale as the reference level, the reference value and parameter estimate are added together. 
For example, if estimates of Linf for two groups are 300 and 5, then adding them gives the Linf of 305 for the 
non-reference group.

\emph{Model Comparisons}

As in the \emph{growthlrt} function based on Kimura (1980), growth curves are tested for differences by using likelihood ratio tests.  
These tests assume homogeneity of variances among groups which is why heterogeneity must be corrected before proceeding. Unlike
the \emph{growthlrt} function, \emph{growthlrt.plus} does not automatically make the comparisons. The user must develop the model structures,
save each oject, and test for differences using the function \emph{compare.lrt.plus}.  Examples are provided below.
}
\value{
  \item{model}{the fitting method and model.}
  \item{results}{list element containing the parameter estimates, standard errors, tests of differences from zero, estimates of
	the maximum likelihood sigma(s), log-likelihood, AIC, BIC, sample sizes, residual degrees of freedom and the residual standard error}  
  \item{variance.covariance}{list element containing the variance covariance matrix.}
  \item{correlation}{list element containing the parameter correlation matrix.}

  \item{residuals}{list element containing the raw and standardized residuals from the model fit.}
  \item{fitted}{list element containing the fitted values from the model fit.}
  \item{convergence}{list element containing convergence information from the \emph{optim} fit.}
  \item{model_comp_df}{list element containing model degrees of freedom used in model comparison.}
  \item{type}{list element containing object type.}

}
\references{Kimura, D. K. 1990. Testing nonlinear regression parameters under heteroscedastic,
 normally-distributed errors. Biometrics 46: 697-708.

Pinheiro, J. C. and D. M. Bates. 2004. Mixed-Effects Models in S and S-PLUS. Springer New York, New York. 528 p.
}

\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}
}
\seealso{\code{\link{growthlrt}} \code{\link{compare.lrt.plus}}}
\examples{
\dontrun{

#### This example produces the same results as the example in 
#### the \emph{growthlrt} helpfile

data(Kimura)

##H0 Model - Different Linfs, Ks and tos for each sex
Ho<-growthlrt.plus(length~Linf*(1-exp(-K*(age-t0))),data=Kimura,
               params=list(Linf~sex,K~sex,t0~sex),
               start=list(Linf=c(60,10),K=c(0.3,0.1),t0=c(0.5,0.05)))

##H1 Model - Same Linfs
H1<-growthlrt.plus(length~Linf*(1-exp(-K*(age-t0))),data=Kimura,
                   params=list(Linf~1,K~sex,t0~sex),
                   start=list(Linf=c(60),K=c(0.3,0.1),t0=c(0.5,0.05)))

##H2 Model - Same Ks
H2<-growthlrt.plus(length~Linf*(1-exp(-K*(age-t0))),data=Kimura,
                   params=list(Linf~sex,K~1,t0~sex),
                   start=list(Linf=c(60,10),K=c(0.3),t0=c(0.5,0.05)))

##H3 Model - Same t0s
H3<-growthlrt.plus(length~Linf*(1-exp(-K*(age-t0))),data=Kimura,
                   params=list(Linf~sex,K~sex,t0~1),
                   start=list(Linf=c(60,10),K=c(0.3,0.1),t0=c(0.5)))

##H4 Model - Same Linf, K and t0
H4<-growthlrt.plus(length~Linf*(1-exp(-K*(age-t0))),data=Kimura,
                   params=list(Linf~1,K~1,t0~1),
                   start=list(Linf=60,K=0.3,t0=0.5))

compare.lrt.plus(Ho,H1)
compare.lrt.plus(Ho,H2)
compare.lrt.plus(Ho,H3)
compare.lrt.plus(Ho,H4)

####This example is Case 2 from Kimura (1990;page 703) and uses the SFR paramterization of the 
#### von Bertalanffy growth equation.

data(AtkaMack)

alt_hypoth_2<-growthlrt.plus(len~lmin+(lmax-lmin)*((1-k^(m-1))/(1-k^(n-1))), 
                   data=AtkaMack,
                   params=list(lmin~sex,lmax~sex,k~sex),
                   within_grp_var=~sex,
                   start=list(lmin=c(26,-2),lmax=c(41,-2),k=c(0.737,0.05)))

null_hypoth_2<-growthlrt.plus(len~lmin+(lmax-lmin)*((1-k^(m-1))/(1-k^(n-1))),
                   data=AtkaMack,
                   params=list(lmin~1,lmax~1,k~1),
                   within_grp_var=~sex,
                   start=list(lmin=c(26),lmax=c(41),k=c(0.737)))

compare.lrt.plus(alt_hypoth_2,null_hypoth_2)

}
}
\keyword{misc}