\name{bonito}
\alias{bonito}
\alias{bonito}

\docType{data}
\title{Data from an age and growth study of the pacific bonito.}
\description{
Growth increment data derived from tagging experiment on Pacific bonito
(Sarda chiliensis) used to illustrate Francis's maximum likelihood
method estimation of growth and growth variability (1988), often refered
to as GROTAG.
}
\usage{data(bonito)}
\format{
  A data frame with 138 observations on the following 4 variables.
  \describe{
    \item{\code{T1}}{a numeric vector describing the release date}
    \item{\code{T2}}{a numeric vector describing the recovery date}
    \item{\code{L1}}{a numeric vector describing the length at release
      in millimeters}
    \item{\code{delta.L}}{a numeric vector describing the length
      increment, in millimeters}
  }
}
\details{
Note that Francis (1988) has discarded 2 records from the original
dataset collected by Campbell et al. (1975).
}

\source{
\describe{\item{1}{Francis, R.I.C.C., 1988. Maximum likelihood estimation of growth and growth variability from tagging data. New Zealand Journal of Marine and Freshwater Research, 22, p.42--51.}
\item{2}{Campbell, G. & Collins, R., 1975. The age and growth of the Pacific bonito, Sarda chiliensis, in the eastern north Pacific. Calif. Dept. Fish Game, 61(4), p.181-200.}
}
}

\keyword{datasets}

