\name{growth_LEP}
\alias{growth_LEP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A flexible maximum likelihood approach for fitting growth curves to tag-recapture data
}
\description{
Estimation of von Bertanffy growth parameters from tag-recapture data (Laslett et al. 2002)
}
\usage{
growth_LEP(l1=NULL,l2=NULL,dt=NULL,measurer = NULL,
           gmodel=1,use_parameter_boundaries=T,graphs=T,
           K_start_bounds=list(K1=NULL,lower_K1=0,upper_K1=Inf,
                 K2=NULL,lower_K2=0,upper_K2=Inf), 
           mu_Linf_start_bounds=list(mu_Linf=NULL,lower_mu_Linf=0,upper_mu_Linf=Inf,
                 sigma_mu_Linf=NULL, lower_sigma_mu_Linf=0,upper_sigma_mu_Linf=Inf),
           A_start_bounds=list(mean_age=NULL,lower_mean_age=0,upper_mean_age=Inf,
                 sigma_age=NULL,lower_sigma_age=0,upper_sigma_age=Inf),
           resid_error_start_bounds=list(sigma_resid=NULL,lower_sigma_resid=0,
                  upper_sigma_resid=Inf),
            measurer_error_start_bounds=list(use_measurer=F,sigma_measure=NULL,
                  lower_sigma_measure=0,upper_sigma_measure=Inf),
           vb_log_k_parms=list(alpha=NULL,lower_alpha=0,upper_alpha=Inf,fix_beta=T,
                 beta=NULL,lower_beta=0,upper_beta=Inf),
           nlminb.control=list(eval.max=10000,iter.max=10000,trace=10),
           tmb.control=list(maxit=10000,trace=FALSE))
}

\arguments{
 \item{l1}{vector of release lengths of tagged fish.}
  \item{l2}{vector of recapture lengths.}
  \item{dt}{vector of time increment between tagging and recapture.}
  \item{measurer}{vector of integers specifying the recapturer type for each row: scientist = 0; fisherperson=1). Not required. Default = NULL.}
  \item{gmodel}{model to fit. 1 = standard von Bertalannfy growth model; 2 = VB log k model of Laslett et al. (2002). Default=1}
  \item{use_parameter_boundaries}{use parameter boundary values (T/F). Applies to all parameters estimated in the model. Default=T.}
  \item{graphs}{plot the observed values of l1 and l2 and the fitted growth curve model versus \emph{Af}(the corrected measures of A; Lasett et al. 2004).  Residuals plots (observed versus fitted) are also provided.}
  \item{K_start_bounds}{list of starting values (K1 and K2), lower(lower_K1 and lower_K2) and upper (upper_K1 and upper_K2) parameter boundaries for K1 and K2. If gmodel = 1, only K1 values are used.}
  \item{mu_Linf_start_bounds}{list of starting, lower and upper boundary values for estimated parameters \emph{mu_Linf} and \emph{sigma_mu_Linf}.}
  \item{A_start_bounds}{list of starting, lower and upper boundaries values for estimated \emph{mean_age}, and starting, lower and upper boundaries values for \emph{sigma_age}, both used  to define the log-normal random effects distribution of A.}
  \item{resid_error_start_bounds}{list of starting, lower and upper boundary values for the estimated residual (measurement) error parameter \emph{sigma_resid}} 
  \item{measurer_error_start_bounds}{list of starting, lower and upper boundary values for the estimated measurer error parameter \emph{sigma_measure}.  Specify \emph{use_measurer}=T to estimate the parameter.  Default is F.}
  \item{vb_log_k_parms}{ If gmodel=2, a list of starting, lower and upper values for estimated parameters \emph{alpha} and \emph{beta}. To fix beta to a constant value, specify \emph{fix_beta}=T and enter a fixed value in \emph{beta}} 
  \item{nlminb.control}{controls for the \emph{nlminb} function. See function \emph{nlminb} for more information.}
  \item{tmb.control}{controls for the \emph{TMB} function. See package \emph{TMB} for more information.}
}

\details{
The von Bertalanffy growth model or the VB log k model of Laslett et al. (2002) is fitted to tag release-capture lengths and times-at-a-large data following Laslett et al. (2002).
The distribution of A is assumed log-normal. 
In addition, adjustments are made to A (age) following Laslett et al. (2004) to correct bias which permits simple graphical checking of the fitted growth curve model. 
If argument \emph{graph} = TRUE, plots of l1 and l2 observed versus predicted, and residuals are created for checking model fit.
Refer to Laslett et al. (2002) for more details.
}
\value{
List containing the parameter_estimates, AIC, random effects A, the original predicted values, the original model residuals, results of the adjustment of A (Af, predicted and residuals for l1 and l2 used for plotting (see Laslett et al., 2004)),
 and convergence statistics (from nlminb; convergence=0 is successful convergence). 
}
\references{

Lasett, G. M., J. P. Eveson and T. Polacheck. 2002. A flexible maximum likelihood approach for fitting growth curves to tag-recapture data.
Canadian Journal of Fisheries and Aquatic Sciences 59: 976-986.

Lasett, G. M., J. P. Eveson and T. Polacheck. 2004. Estimating the age at capture in capture-recapture studies of fish growth. Australian and New Zealand Journal of Statistics 46: 59-66.

}
\author{
Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}

}
\note{
Paige Eveson of CSIRO Marine Research kindly provided the R code for calculating \emph{Af} based on Lasett et al., 2004.
}

\seealso{
\code{\link{growhamp}} \code{\link{grotag}} \code{\link{grotagplus}}
}
\examples{
 \dontrun{
	data(lepdata)
  growth_LEP(l1=lepdata$l1,l2=lepdata$l2,dt=lepdata$dt,measurer=NULL,
               gmodel=1,use_parameter_boundaries=T,graphs=T,
               K_start_bounds=list(K1=0.2,lower_K1=0,upper_K1=Inf,K2=0.12,lower_K2=0,
                  upper_K2=Inf), 
               mu_Linf_start_bounds=list(mu_Linf=189.624,lower_mu_Linf=0,upper_mu_Linf=Inf,
                  sigma_mu_Linf=11.032,lower_sigma_mu_Linf=0, upper_sigma_mu_Linf=Inf),
               A_start_bounds=list(mean_age=1.76,lower_mean_age=0,upper_mean_age=Inf,
                   sigma_age=0.165,lower_sigma_age=0,upper_sigma_age=Inf),
               resid_error_start_bounds=list(sigma_resid=3.547,lower_sigma_resid=0,
                    upper_sigma_resid=Inf),
               measurer_error_start_bounds=list(use_measurer=F,sigma_measure=3.547,
                     lower_sigma_measure=0,upper_sigma_measure=Inf),
               vb_log_k_parms=list(alpha=2.955,lower_alpha=0,upper_alpha=Inf,fix_beta=T,
                      beta=30,lower_beta=0,upper_beta=30),
                      nlminb.control=list(eval.max=10000,iter.max=10000,trace=10),
                      tmb.control=list(maxit=10000,trace=FALSE))
 }
}
\keyword{misc}
