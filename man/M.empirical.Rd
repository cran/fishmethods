\name{M.empirical}
\alias{M.empirical}
\title{Estimation of Natural Mortality Rates from Life History Parameters}
\description{
 The approaches of Pauly (1980), Hoenig (1983), Alverson and Carney (1975), Roff (1984), Gunderson and Dygert (1988),
Petersen and Wroblewski (1984), Lorenzen (1996), Gislason et al. (2010), Then et al. (2015), Brey (1999) and 
 Charnov et al. (2013) are encoded for estimation of natural mortality (M).
}
\usage{
M.empirical(Linf = NULL, Winf = NULL, Kl = NULL, Kw = NULL,
 TC = NULL, tmax = NULL, tm = NULL, GSI = NULL, Wdry = NULL,
 Wwet = NULL, Bl = NULL, TK = NULL, BM = NULL, L = NULL, method = c(1, 2, 
3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
}
\arguments{
  \item{Linf}{Length-infinity value from a von Bertalanffy growth curve (total length-cm).}
  \item{Winf}{Weight-infinity value from a von Bertalanffy growth curve (wet weight-grams).}
  \item{Kl}{Kl is the growth coefficient (per year) from a von Bertalanffy growth curve for length.}
  \item{Kw}{Kw is the growth coefficient (per year) from a von Bertalanffy growth curve for weight.}
  \item{TC}{the mean water temperature (Celsius) experienced by the stock.}
  \item{tmax}{the oldest age observed for the species.}
  \item{tm}{the age at maturity.}
  \item{GSI}{gonadosomatic index (wet ovary weight over wet body weight).}
  \item{Wdry}{total dry weight in grams.}
  \item{Wwet}{total wet weight at mean length in grams.}
  \item{Bl}{body length in cm.}
  \item{TK}{mean temperature (Kelvin).}
  \item{BM}{maximum body mass (kJ - kiloJoules)}
  \item{L}{fish length along the growth trajectory}
  \item{method}{vector of method code(s).  Any combination of methods can employed. \code{1}= Pauly (1980)
 length equation - requires Linf, Kl, and TC; \code{2}= Pauly (1980) weight equation - requires Winf, Kw, and TC; 
\code{3}= Hoenig (1983) joint equation - requires tmax; \code{4}= Alverson and Carney (1975) - requires Kl and tmax; 
\code{5}= Roff (1984) - requires Kl and tm; \code{6}= Gunderson and Dygert (1988) - requires GSI; \code{7}= Peterson
 and Wroblewski (1984) - requires Wdry; \code{8}= Lorenzen (1996) - requires Wwet;
 \code{9}= Gislason et al. (2010) - requires Linf, K and Bl; \code{10}= Then et al. (2015) tmax - requires tmax;
\code{11}= Then et al. (2015) growth  - requires Kl and Linf.
\code{12}= Brey (1999) - requires tmax, TK, and BM.
\code{13}= Charnov et al (2013) - requires Linf, Kl, and L.

  }
}
\details{
 Please read the references below for details about equations. Some estimates of M will not be valid for
certain fish groups.
}
\value{A matrix of M estimates.
}
\references{Alverson, D. L. and M. J. Carney. 1975. A graphic review of the growth and decay of population cohorts. J. Cons. Int. Explor. Mer 36: 133-143.

Brey, T. 1999. Growth performance and mortality in aquatic macrobenthic invertebrates. Advances in Marine Biology 35: 155-223.

Charnov, E. L., H. Gislason, J. G. Pope. 2013. Evolutionary assembly rules for fish life histories. Fish and Fisheries 14: 213-224.

Gislason, H., N. Daan, J. C. Rice, and J. G. Pope. 2010. Size, growth, temperature and the natural mortality of marine fish. Fish and Fisheries 11: 149-158.

Gunderson, D. R. and P. H. Dygert. 1988. Reproductive effort as a predictor of natural mortality rate. J. Cons. Int. Explor. Mer 44: 200-209.

Hoenig, J. M. 1983. Empirical use of longevity data to estimate mortality rates. Fish. Bull. 82: 898-903.

Lorenzen, K. 1996. The relationship between body weight and natural mortality in juvenile and adult fish: a comparison of natural ecosystems and aquaculture. J. Fish. Biol. 49: 627-647.

Pauly, D. 1980. On the interrelationships between natural mortality, growth parameters, and mean environmental temperature in 175 fish stocks. J. Cons. Int. Explor. Mer: 175-192.

Peterson, I. and J. S. Wroblewski. 1984. Mortality rate of fishes in the pelagic ecosystem. Can. J. Fish. Aquat. Sci. 41: 1117-1120.

Roff, D. A. 1984.  The evolution of life history parameters in teleosts. Can. J. Fish. Aquat. Sci. 41: 989-1000. 

Then, A. Y., J. M. Hoenig, N. G. Hall, D. A. Hewitt. 2015. Evaluating the predictive performance of empirical estimators of natural mortality rate using information on over 200 fish species. 
ICES J. Mar. Sci. 72: 82-92.

 }
\note{Original functions for the Pauly (1980) length equation and the Hoenig (1983) fish equation were provided by Michael H. Prager, National Marine Fisheries Service, Beaufort, North Carolina.}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@mass.gov}}
\examples{
M.empirical(Linf=30.1,Kl=0.31,TC=24,method=c(1))
}
\keyword{misc}