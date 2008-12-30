\name{catchseries}
\alias{catchseries}
\title{Summary of Catch Estimates from the Marine Recreational Fisheries Statistics Survey}
\description{
  Species-specific estimates of harvest weight and numbers of fish, released fish, and total catch are summarized for any combination of year(s), state(s),
 bimonthly wave(s), fishing mode(s), and fishing area(s) from the Marine Recreational Fisheries Statistics Survey (MRFSS) data.
}

\usage{
catchseries(estdir = NULL, species = NULL, state = NULL,
 byst = 1, wave = NULL, bywave = 0, mode = NULL, 
 bymode = 0, area = NULL, byarea = 0, styr = NULL,
 endyr = NULL)
}
\arguments{
  \item{estdir}{the path and main directory (in quotes) under which catch/effort sub-directories are stored.}
  \item{species}{10-digit NODC species code. A complete list of species codes is provided at the MRFSS website}
  \item{state}{state(s) code of data to include.  A complete list of state codes is provided at the MRFSS website.
Atlantic and Gulf States: 23= ME, 33= NH, 25= MA, 44= RI, 9= CT, 
      36= NY, 34= NJ, 10= DE, 24= MD, 51= VA,37= NC, 45= SC, 
      13= GA, 12= FL (121= East FL, 122= West FL), 1= AL,
      28=M S, and 22= LA. Any combination of states can be included in c().
      North Atlantic= c(9, 23, 25, 33, 44), 
	Mid Atlantic= c(10, 24, 34, 36, 51), 
	South Atlantic= c(13, 37, 45, 121), 
	Gulf of Mexico= c(1, 22, 28, 122).
 }

 \item{byst}{numeric value. 0 = combine data from states, 1 = show data on state-by-state basis. Default=1 }
 \item{wave}{bimothly wave(s) of data to include. Bimothly Waves:
	1= Jan-Feb, 2= Mar-Apr, 3= May-Jun, 4= Jul-Aug, 5= Sept-Oct, 6= Nov-Dec.
     Any combination of waves can be included in c(). 
 }


  \item{bywave}{numeric value. 0 = combine data across waves, 1 = show data on wave-by-wave basis. Default=0 }
\item{mode}{ the fishing mode(s) of data to include. Valid codes will depend on year.
Years 1982-2004: 3= Shore, 6= Party & Charter Boats, and 7= Private & Rental Boats.
 Years 2005-present: 3= Shore, 4= Party Boat, 5= Charter Boat, and 7= Private & Rental Boats.           
Any combination of modes can be included in c().          
}
 \item{bymode}{numeric value. 0 = combine data across modes, 1 = show data on mode-by-mode basis. Default=0 }
 \item{area}{area code of data to include. Area codes: 1= State Territorial Seas (Ocean<=3 mi excluding Inland), 
2= Exclusive Economic Zone (Ocean>3 mi), 
3= Ocean <=10 mi West FL and TX, 4= Ocean > 10 mi West FL and TX, 5= Inland, 
6= Unknown. Any combination of areas can be included in c().
}

  \item{byarea}{numeric value. 0 = combine data across areas, 1 = show data on area-by-area basis. Default=0 }
  \item{styr}{the starting year of data to include.}
  \item{endyr}{the ending year of data to include.}
}
\details{
 Catch/effort estimates must be extracted using function \code{extractMRFSS} before summarization. Harvest, release and weight estimates are 
always summed by year regardless of whether other factors (state, mode, wave, and area) are aggregated or not. 
Additional statistics (standard error and proportional standard error) are provided for each estimate.
Definition of estimates can be found at the MRFSS wesbite.
}
\value{
  \item{Details}{a list element summarizing the selected arguments.}
  \item{Results}{a list element containing the table of summary statistics (i.e., estimate, standard errors, and proportional
  standard errors) of harvest and release numbers and harvest weight (kg).}
}
\references{Marine Recreational Fisheries Statistics Survey

\url{http://www.st.nmfs.noaa.gov/st1/recreational/index.html}
}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\seealso{\code{\link{extractMRFSS}}}
\examples{
\dontrun{
outpt<-catchseries(estdir="C:/Temp",species=8835250101,
  state=c(37,45,13,121),mode=c(4,5,6),bymode=1,wave=c(1:6),
  area=c(1,2,3,4,5,7), styr=2000,endyr=2007)}
}
\keyword{ misc}