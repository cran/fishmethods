\name{catchpertrip}
\alias{catchpertrip}
\title{Catch-Per-Trip Analysis for Marine Recreational Fisheries Statistics Survey Data}
\description{
  Catch-per-trip analysis for Marine Recreational Fisheries Statistics Survey data following the Chapter 4 of Recreational Fisheries Data User's Manual. 
}
\usage{
catchpertrip(intdir = NULL, estdir = NULL, species = NULL,
 state = NULL, mode = NULL, wave = NULL, styr = NULL,
 endyr = NULL)
}
\arguments{
  \item{intdir}{the path and main directory (in quotes) under which raw intercept sub-directories are stored.}
  \item{estdir}{the path and main directory (in quotes) under which catch/effort sub-directories are stored.}
  \item{species}{10-digit NODC species code. A complete list of species codes is provided at the MRFSS website.}
 \item{state}{state(s) code of data to include.  A complete list of state codes is provided at the MRFSS website.
Atlantic and Gulf States: 23= ME, 33= NH, 25= MA, 44= RI, 9= CT, 
      36= NY, 34= NJ, 10= DE, 24= MD, 51= VA, 37= NC, 45= SC, 
      13= GA, 12= FL (121= East FL, 122= West FL), 1= AL,
      28= MS, and 22= LA. Any combination of states can be included in c().
      North Atlantic= c(9, 23, 25, 33, 44), 
	Mid Atlantic= c(10, 24, 34, 36, 51), 
	South Atlantic= c(13, 37, 45, 121), 
	Gulf of Mexico= c(1, 22, 28, 122).
 }
  \item{wave}{bimothly wave(s) of data to include. Bimothly Waves:
	1= Jan-Feb, 2= Mar-Apr, 3= May-Jun, 4=J ul-Aug, 5= Sept-Oct, and 6= Nov-Dec.
     Any combination of waves can be included in c(). 
 }
  \item{mode}{ the fishing mode(s) of data to include. Valid codes will depend on year.
Years 1982-2004: 3= Shore, 6= Party & Charter Boats, 7= Private & Rental Boats.
 Years 2005-present: 3= Shore, 4= Party Boat, 5= Charter Boat, and 7= Private & Rental Boats.           
Any combination of modes can be included in c().
           
}
  \item{styr}{ the starting year of data to include.}
  \item{endyr}{  the ending year of data to include.}

}
\details{
 Raw intercept data and catch/effort estimates must be extracted using function \code{extractMRFSS} 
before catch-per-trip analysis can be conducted. When performing catch-per-trip analyses within a state/mode/wave 
stratum, reweighting of data is not necessary since intercept sampling is random within each stratum and 
an assumption can be made that the true effort distributions are represented. However, if catch-per-trip analyses are
performed among state/mode/wave strata, data must be reweighted prior to pooling among strata.  This function
properly reweights and combines data for catch-per-trip analyses.
}
\value{
   \item{Results}{Dataframe containing catch level (CATCH), estimated number of trips (TRIPS),
 intercept frequency (CF) and proportion of successful
 trips(PROP.TRIP)}
  }
\references{Marine Recreational Fisheries Statistics Survey - Recreational Fisheries Data User's Manual 

\url{http://www.st.nmfs.noaa.gov/st1/recreational/SurveyMaterials.html}}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{extractMRFSS}}}
\examples{
## Example of a typical specification - for Mid-Atlantic
## black sea bass wave 4, all modes in 1992. Figure 11 of Chapter 4 
## in User Manual
\dontrun{catchpertrip(intdir="C:/Temp",estdir="C:/Temp",species=8835020301,
state=c(10,24,34,36,51),mode=c(3,6,7),wave=c(4), styr=1992,endyr=1993)}
}
\keyword{misc}
