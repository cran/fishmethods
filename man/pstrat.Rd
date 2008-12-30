\name{pstrat}
\alias{pstrat}
\title{Post-stratification of Marine Recreational Fisheries Statistics Survey Estimates of Catch}
\description{
  Post-stratified catch and variance estimates are generated from the original intercept data and catch & effort estimates provided by 
the Marine Recreational Fisheries Statistics Survey (MRFSS).
}
\usage{pstrat(intdir = NULL, estdir = NULL, pstdir = NULL, state = NULL,
 year = NULL, stwave = NULL, endwave = NULL, psfactor = NULL)
}
\arguments{
  \item{intdir}{the path and main directory (in quotes) under which raw intercept sub-directories are stored.}
  \item{estdir}{the path and main directory (in quotes) under which catch/effort sub-directories are stored.}
  \item{pstdir}{the path and main directory (in quotes) under which post-stratified estimates will be stored. The program will
           store the post-stratification estimates as a .csv file under a subdirectory named \emph{pest}\code{year}.}
 \item{year}{year of data to include.} 
   \item{state}{state(s) code of data to include.  A complete list of state codes is provided at the MRFSS website.
Atlantic and Gulf States: 23= ME, 33= NH, 25= MA, 44= RI, 9= CT, 
      36= NY, 34= NJ, 10= DE, 24= MD, 51= VA, 37= NC, 45= SC, 
      13= GA, 12= FL (121= East FL, 122= West FL), 1=AL,
      28= MS, and 22= LA. Any combination of states can be included in c().
      North Atlantic= c(9, 23, 25, 33, 44), 
	Mid Atlantic= c(10, 24, 34, 36, 51), 
	South Atlantic= c(13, 37, 45, 121), 
	Gulf of Mexico= c(1, 22, 28, 122).
 }
  \item{stwave}{bimothly wave(s) of data to include. Bimothly Waves:
	1= Jan-Feb, 2= Mar-Apr, 3= May-Jun, 4= Jul-Aug, 5= Sept-Oct, and 6= Nov-Dec.
     Any combination of waves can be included in c(). 
 }
 
  \item{endwave}{last bimonthly wave of data to include.}
  \item{psfactor}{name of dataframe with poststratification variable(s)}
}
\details{
   Raw intercept data and catch/effort estimates must be extracted using function \code{extractMRFSS}.  Post-stratification methodology follows 
the SAS post-stratification program available from MRFSS. A dataframe must be created that contains
the MRFSS variable(s) linked to the variable AREA.G which identifies the post-stratification scheme.  

As a first example, the State of Massachusetts post-stratifies catch estimates by county to examine regional fishing 
harvest.  The MRFSS variable CNTY is used to group the regions.  To code the post-stratification variable, a dataframe 
containing CNTY numbers and the associated grouping labels under variable \code{AREA.G} is created as follows:

\preformatted{post<-data.frame(CNTY=c(9,21,25,23,1,7,5,18,17,27),
AREA.G=c("A","B","C","D","E","F","G","H","Y","Y"))
}

Each \code{AREA.G} label is matched to a specific county.  \code{AREA.G} must be the name of the grouping label. More than one variable (other than CNTY) can be used to parse
the data into more specific strata by creating an additional column in the dataframe.


As a second example, the State of Florida post-stratifies catch estimates into four regions (Northeast, Southeast,
Southwest and Northwest Florida) related to spotted seatrout (Cynoscion nebulosus) management. The coding for 
the four regions uses CNTY to create groupings:

\preformatted{post<-data.frame(CNTY=c(89,31,109,35,19,107,127,9,61,111,85,
                        99,11,25,87,21,51,71,15,27,115,81,
                        57,103,101,53,17,75,1,29,123,65,129,37,
                        77,45,5,133,131,91,113,33
                        ),
                       AREA.G=c("N","N","N","N","N","N",
                       "S","S","S","S","S","S","S",
                       "S","S","S","S","S","S","S",
                       "S","S","S","S","N","N","N",
                       "N","N","N","N","N","N","N",
                       "N","N","N","N","N","N","N",
                       "N")) 
}      

Notice that there are not four \code{AREA.G} labels for the different areas because post-stratification is automatically performed 
over the sub-region variable (Florida East Coast SUB.REG=6 and Florida West Coast is SUB.REG=7).  For states with only one SUB.REG code,
four \code{AREA.G}labels would need to be created. 

The AREA.S variable, which is the AREA.G and AREA.X labels combined, is created in the output file and designates the
post-stratification levels.

All variables coded found in the MRFSS Recreational Fisheries Data User's Manual.

}
\value{
  Post-stratified estimates from each year and wave are written to comma-delimited (.csv) files and are stored under
the \code{pestdir} directory as \code{AG_}\emph{yyyywv}\code{.csv} (e.g., \code{AG_}\code{20075.csv}). 
}
\references{Marine Recreational Fisheries Statistics Survey - Recreational Fisheries Data User's Manual

 \url{http://www.st.nmfs.noaa.gov/st1/recreational/SurveyMaterials.html}}
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\note{The MRFSS post-stratification program written in SAS by Gerry Gray et al. of the National Marine Fisheries Service was converted to R-code.
}
\seealso{\code{\link{extractMRFSS}}}
\examples{
## An example of code - not a working example
\dontrun{
post<-data.frame(CNTY=c(9,21,25,23,1,7,5,18,17,27),
            AREA.G=c("A","B","C","D","E","F","G","H","Y","Y"))

pstrat(intdir="C:/Temp", estdir="C:/Temp",pstdir="C:/Temp",state=25,
  year=2007,stwave=2,endwave=6,psfactor=post)
}
}
\keyword{misc}