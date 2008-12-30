\name{lengthfreq}
\alias{lengthfreq}
\title{Calculate Length Frequencies from Marine Recreational Fisheries Statistics Survey Data}
\description{
 A combined length frequency table for a species is generated for any combination of year(s), state(s), bimonthly wave(s),
 fishing mode(s), and fishing area(s) from the Marine Recreational Fisheries Statistics Survey (MRFSS) data.  The analytical steps 
 given in Recreational Fisheries Data User's Manual are followed.
}
\usage{
lengthfreq(intdir = NULL, estdir = NULL, species = NULL,
 state = NULL, wave = NULL, mode = NULL, area = NULL,
 styr = NULL, endyr = NULL, conveq = FALSE,
 parms = c(0, 1))
}
\arguments{
  \item{intdir}{the path and main directory (in quotes) under which raw intercept sub-directories are stored.}
  \item{estdir}{the path and main directory (in quotes) under which catch/effort sub-directories are stored.}
  \item{species}{10-digit NODC species code. A complete list of species codes is provided at the MRFSS website}
  \item{state}{state(s) code of data to include.  A complete list of state codes is provided at the MRFSS website.
Atlantic and Gulf States: 23= ME, 33= NH, 25= MA, 44= RI, 9= CT, 
      36= NY, 34= NJ, 10= DE, 24= MD, 51= VA,37= NC, 45= SC, 
      13= GA, 12= FL (121= East FL, 122= West FL), 1=AL,
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
Years 1982-2004: 3= Shore, 6= Party & Charter Boats, and 7= Private & Rental Boats. Years 2005-present:
 3= Shore, 4= Party Boat, 5= Charter Boat, and 7= Private & Rental Boats.           
Any combination of modes can be included in c().
           
}

 \item{area}{area code of data to include. Valid codes are: 1= State Territorial Seas (Ocean<= 3 mi excluding Inland),
 2= Federal Exclusive Economic Zone (Ocean > 3 mi),
3= Ocean <=10 mi West FL and TX, 4= Ocean > 10 mi West FL and TX, 5= Inland, and 6= Unknown. Any combination of areas
 can be included in c().
}

  \item{conveq}{logical to indicate whether a length conversion equation will be provided. If TRUE, provide the intercept 
and slope in \code{parms} argument below. Conversion equation
is used to convert MRFSS fork length (mm) to other length measurements.}
  \item{parms}{vector containing the intercept (position 1) and slope  (position 2) of the conversion equation}
  \item{styr}{ the starting year of data to include.}
  \item{endyr}{  the ending year of data to include.}
}
\details{
 Raw intercept data and catch/effort estimates must be extracted using function \code{extractMRFSS} before length frequencies can be created.
 Length frequencies are produced with one-inch length groups identical to the MRFSS website.  Before creation of the length frquency table,
 fork length (in millimeters) is converted to inches. The intercept and slope parameters of a conversion equation can be introduced to convert fork length into 
 other length measurements before creating the length frequencies. Length data and harvest estimates from the catch/effort files are extracted using the argument values. 
 Length frequencies are created for each state/mode/wave/area fished stratum and are weighted with matching harvest estimates before
 being combined across strata.
}
\value{
  \item{Details}{a list element summarizing the selected arguments.}
  \item{Results}{a list element containing the length frequency table with one-inch length group, numbers-at-length, and proportions-at-length}
}
\references{Marine Recreational Fisheries Statistics Survey - Recreational Fisheries Data User's Manual

 \url{http://www.st.nmfs.noaa.gov/st1/recreational/SurveyMaterials.html}
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{extractMRFSS}}}
\examples{
##Example of a typical specification - not a working example
\dontrun{
dodo<-lengthfreq(intdir="C:/Temp",estdir="C:/Temp",species=8835250101,
 state=c(25),mode=c(3,4,5,7),wave=c(1,2,3,4,5,6),area=c(1,2,3,4,5,7),
 styr=2007,endyr=2007)}
}
\keyword{misc}
