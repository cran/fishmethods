\name{lifetable}
\alias{lifetable}
\title{Life Table Construction}
\description{
  Life tables are constructed from either numbers of individuals of a cohort alive at the
start of an age interval (nx) or number of individuals of a cohort dying during the age interval (dx).

}
\usage{
lifetable(age = NULL, numbers = NULL, r = NULL, type = 1)
}
\arguments{
  \item{age}{ vector of age intervals (e.g., 0 to maximum cohort age).}
  \item{numbers}{number of individual alive (nx) or dead (dx) }
  \item{r}{known rate of increase (r) for methods 3 and 4}
  \item{type}{numeric value of method to use to calculate life table.

 1 = Age at death recorded directly and no assumption made about population stability or stability of age structure - Method 1 in Krebs (1989).
 2 = Cohort size recorded directly and and no assumption made about population stability or stability of age structure - Method 2 in Krebs (1989).
 3 = Ages at death recorded for a population with stable age distribution and known rate of increase - Method 5 in Krebs (1989).
 4 = Age distribution recorded for a population with a stable age distribution and known rate of increase - Method 6 in Krebs (1989).
}
}
\details{
 Following Krebs (1989:413-420),  standard life tables are calculated given age intervals and either cohort size or deaths.
X=age interval, nx=number of individuals of a cohort alive at the start of age interval X, lx = proportion of individuals surviving at the start of
age interval X, dx = number of individuals of a cohort dying during the age interval X,
qx=finite rate of mortality during the age interval X to X+1, px=finite rate of survival during the age interval X to X+1,
ex=mean expectation of life for individuals alive at start of age X. 
For method 5, dx is corrected for population growth by \emph{dx'=dx*exp(r*x)} and in method 6, nx is corrected for the same by \emph{nx*e(r*x)}.
See Krebs for formulae.
}
\value{Dataframe containing life table values.
}
\references{Krebs, C. J. 1989. \emph{Ecological Methodologies}. Harper and Row, New York, NY. 654 p.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
data(buffalo)
lifetable(age=buffalo$age,numbers=buffalo$nx,type=2)
}
\keyword{misc}
