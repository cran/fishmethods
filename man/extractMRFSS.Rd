\name{extractMRFSS}
\alias{extractMRFSS}
\title{Extraction of Marine Recreational Fisheries Statistics Survey Data from SAS Transport Files}
\description{MRFSS intercept and catch/effort data in SAS transport files (.xpt) are extracted using 
  package \code{foreign} and are saved as .csv files under data and year-specific subdirectories.
}
\usage{
extractMRFSS(indir = NULL, outdir = NULL, type = NULL,
 state = NULL, styr = NULL, endyr = NULL)
}
\arguments{
  \item{indir}{the path and directory (in quotes) under which the SAS transport files are stored.}
  \item{outdir}{the path and main directory (in quotes) under which data- and year-specific subdirectories will be created
               and .csv files will be stored.}
  \item{type}{the type of data to extract. 1 = raw intercept data, 2= catch/effort estimates.}
  \item{state}{the state code(s) designating which state(s) data to extract.  If multiple states, separate state codes with commas within the 
 concatentation function.  If data from all states are desired, do not include the \code{state=} argument. State codes are 
 found on the MRFSS website.}
  \item{styr}{the starting year of data to extract.}
  \item{endyr}{the ending year of data to extract.}
}
\details{Download the zipped intercept data and catch/effort estimates from the MRFSS website \url{http://www.st.nmfs.noaa.gov/st1/recreational/index.html}. Unzip/extract the 
SAS transport files (.xpt) to the \code{indir} directory. Using this function, specify the \code{indir} and \code{outdir} directories, which data type to extract, the state code(s), 
and the starting and ending years (from .xpt files) of data to be extracted. If the \code{outdir} directory is not specified, the files will be stored under the 
\code{indir} directory.
}
\value{
  Raw intercept data and catch & effort estimates saved as .csv files.
 }
\references{\url{http://www.st.nmfs.noaa.gov/st1/recreational/index.html}. }
\author{Gary A. Nelson, Massachusetts Division of Marine Fisheries \email{gary.nelson@state.ma.us}}
\examples{
## This is a typical specification, not a working example
\dontrun{extractMRFSS(indir="C:/Temp",outdir="C:/MRFSS",type=1,state=25,
styr=1981,endyr=2007)}
}
\keyword{misc}