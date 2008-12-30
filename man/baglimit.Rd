\name{baglimit}
\alias{baglimit}
\title{Bag Limit Analysis of Marine Recreational Fisheries Statistics Survey Data}
\description{Analysis of Marine Recreational Fisheries Statistics data to determine the effect a simulated bag limit would have on the overall reduction in harvest.
}
\usage{
baglimit(intdir = NULL, estdir = NULL, species = NULL, 
  state = NULL, mode = NULL, wave = NULL, styr = NULL, 
  endyr = NULL, bag = NULL)
}
\arguments{
  \item{intdir}{the path and main directory (in quotes) under which raw intercept sub-directories are stored.}
  \item{estdir}{the path and main directory (in quotes) under which catch/effort sub-directories are stored.}
  \item{species}{10-digit NODC species code. A complete list of species codes is provided at the MRFSS website}
  \item{state}{state(s) code of data to include.  A complete list of state codes is provided at the MRFSS website.
Atlantic and Gulf States: 23= ME, 33= NH, 25= MA, 44= RI, 9= CT, 
      36= NY, 34= NJ, 10= DE, 24= MD, 51= VA,37= NC, 45= SC, 
      13= GA, 12= FL (121=East FL, 122=West FL), 1=AL,
      28= MS, and 22= LA. Any combination of states can be included in c().
      North Atlantic= c(9, 23, 25, 33, 44), 
	Mid Atlantic= c(10, 24, 34, 36, 51), 
	South Atlantic= c(13, 37, 45, 121), 
	Gulf of Mexico= c(1, 22, 28, 122).
 }
  \item{wave}{bimothly wave(s) of data to include. Bimothly Waves:
	1= Jan-Feb, 2= Mar-Apr, 3= May-Jun, 4= Jul-Aug, 5= Sept-Oct, and 6= Nov-Dec.
     Any combination of waves can be included in c(). 
 }
  \item{mode}{ the fishing mode(s) of data to include. Valid codes will depend on year.
Years 1982-2004: 3= Shore, 6= Party & Charter Boats, and 7= Private & Rental Boats.
 Years 2005-present: 3= Shore, 4= Party Boat, 5= Charter Boat, and 7= Private & Rental Boats.           
Any combination of modes can be included in c().
           
}
  \item{styr}{ the starting year of data to include.}
  \item{endyr}{the ending year of data to include.}
  \item{bag}{the bag limit(s) to simulate.  If multiple bag limits, separate values by commas in c().}
}
\details{
Raw intercept and catch & effort estimates must be extracted using function \code{extractMRFSS} before bag limit simulations can be done.
Species-specific intercept data and associated effort estimates for each state, year, wave and mode are extracted.  Harvest numbers with multiple 
contributors (CNTRBTRS) are split evenly among the number of contributors (fractional harvest frequencies may result). The analytical steps given 
in Chapter 7 of \code{Recreational Fisheries Data Users Manual} are followed. Estimated total harvest from the \code{No Bag Limit} table should be close 
to the MRFSS harvest estimate.}
\value{
  \item{No Bag Limit}{list element containing the table of catch, trip,intercept and harvest information with no bag limit.}
  \item{Bag Limit i}{list elements containing the table of catch, trip,intercept and harvest information with each bag limit applied.
        The number of elements will depend on number of bag limit simulations.}
  \item{Results}{the percent reduction in harvest associated with each bag limit.}
  \item{Details}{a list element summarizing the selected arguments.}
}
\references{Marine Recreational Fisheries Statistics Survey - Recreational Fisheries Data User's Manual

 \url{http://www.st.nmfs.noaa.gov/st1/recreational/SurveyMaterials.html}
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{extractMRFSS}}}
\examples{
## Example of a typical specification. Example is for bluefish from New  
## York in 1990 (all waves and modes)- not a working example
\dontrun{
baglimit(intdir="C:/Temp",estdir="C:/Temp",species=8835250101,
 state=c(36),mode=c(3,6,7),wave=c(1,2,3,4,5,6),bag=c(5,10,15),
 styr=1990,endyr=1990)}
}
\keyword{misc}