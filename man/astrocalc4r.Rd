\name{astrocalc4r}
\alias{astrocalc4r}
\title{Solar zenith angles for biological research
}
\description{
This function calculates the solar zenith, azimuth and declination angles, time at sunrise, local noon and sunset, day length, 
and PAR (photosynthetically available radiation, 400-700 nm) under clear
skies and average atmospheric conditions
(marine or continental) anywhere on the surface of the earth based on date, time, and location.}
\usage{
astrocalc4r(day, month, year, hour, timezone, lat, lon, withinput = FALSE, 
seaorland = "maritime", acknowledgment = FALSE)
}
\arguments{
  \item{day}{day of month in the local time zone (integers). Value is
    required. Multiple observations should be
    enclosed with the c() function.}
  \item{month}{month of year in the local time zone (integers). Value is
    required. Multiple observations should be 
    enclosed with the c() function.}
  \item{year}{year in the local time zone (integers).Value is
    required. Multiple observations should be
    enclosed with the c() function.}
  \item{hour}{local time for each observation (decimal hours, e.g. 11:30 PM is 23.5, real numbers). Value is required. 
    Multiple observations should be enclosed with the c() function.}
  \item{timezone}{local time zone in +/- hours relative to GMT to link
    local time and GMT. For example, the difference
    between Eastern Standard Time and GMT is -5 hours. Value is
    required. Multiple observations should be enclosed 
    with the c() function. timezone should include any necessary adjustments for daylight savings time.}
  \item{lat}{Latitude in decimal degrees (0o to 90 o in the northern  hemisphere and -90 o to 0 o degrees in the 
    southern hemisphere, real numbers). For example, 42o 30' N is 42.5 o and 42o 30' S is -42.5o. 
   Value is required. Multiple observations should be enclosed with the c() function.}
  \item{lon}{Longitude in decimal degrees (-0 o to 180 o in the western
    hemisphere and 0o to 180 o in the eastern hemisphere, real numbers). 
     For example, 110o 15' W is -110.25 o and 110o 15' E is 110.25o. Value is required. 
     Multiple observations should be enclosed with the c() function.}
  \item{withinput}{logical:TRUE to return results in a dataframe with
    the input data; otherwise FALSE returns a dataframe with just results. Default is FALSE.}
  \item{seaorland}{text: "maritime" for typical maritime conditions or "continental" for typical continental conditions.
    Users must select one option or the other based on proximity to the ocean or other factors.}
  \item{acknowledgment}{logical: use TRUE to output acknowledgement. Default is FALSE.}
}
\details{
Astronomical definitions are based on definitions in Meeus (2009) and Seidelmann (2006).
The solar zenith angle is measured between a line drawn "straight up" from the center of the earth through the 
observer and a line drawn from the observer to the center of the solar disk. 
The zenith angle reaches its lowest daily value at local noon when the sun is highest. It reaches its maximum value at 
night after the sun drops below the horizon.  The zenith angle and all of the solar variables calculated by
\code{astrocalc4r} depend on latitude, longitude, date and time of day. For example, solar zenith angles measured 
at the same time of day and two different locations would differ due to
differences in location. Zenith angles at the same location and 
two different dates or times of day also differ.

Local noon is the time of day when the sun reaches its maximum elevation and minimum solar zenith angle at the observers location. This 
angle occurs when the leading edge of the sun first appears above, or the trailing edge disappears below 
the horizon (0.83o accounts for the radius of the sun when seen from the earth and for refraction by the atmosphere). 
Day length is the time in hours between sunrise and sunset. Solar
declination and azimuth angles describe the exact position of the sun in
the sky relative to an observer based on an equatorial coordinate system (Meeus 2009). Solar declination 
is the angular displacement of the sun above the equatorial plane. The
equation of time accounts for the relative
position of the observer within the time zone and is provided because it
is complicated to calculate. PAR isirradiance in lux (lx, approximately W m-2) at the surface of the 
earth under clear skies calculated based on the solar zenith angle and assumptions about marine or terrestrial atmospheric 
properties. \code{astrocalc4r} calculates PAR for wavelengths between 400-700 nm. Calculations for other wavelengths
 can be carried out by modifying the code to use parameters from Frouin et al. (1989). Following Frouin et al. (1989), 
PAR is assumed to be zero at solar zenith angles >= 90o although some sunlight may be visible in the sky 
when the solar zenith angle is < 108o. Angles in \code{astrocalc4r} output are in degrees although radians are used 
internally for calculations. Time data and results are in decimal hours (e.g. 11:30 pm = 23.5 h) local time but internal
 calculations are in Greenwich Mean Time (GMT). The user must specify the local time zone in terms of +/- hours relative to GMT to link 
local time and GMT. For example, the difference between Eastern Standard Time and GMT is -5 hours.
The user must ensure that any adjustments for daylight savings time are included in the timezone value. For example, 
timezone=-6 for Eastern daylight time.
}
\value{Time of solar noon, sunrise and sunset, angles of azimuth and zenith, eqtime, declination of sun,
 daylight length (hours) and PAR.
}
\references{
Frouin, R., Lingner, D., Gautier, C., Baker, K. and Smith, R. 1989. A simple analytical formula 
to compute total and photosynthetically available solar irradiance at the ocean surface under 
clear skies. J. Geophys. Res. 94: 9731-9742. 

L. D. Jacobson, L. C. Hendrickson, and J. Tang. 2015. Solar zenith angles for biological research and an expected
catch model for diel vertical migration patterns that affect stock size estimates for longfin inshore squid
 (Doryteuthis pealeii). Canadian Journal of Fisheries and Aquatic Sciences 72: 1329-1338.

Meeus, J. 2009. Astronomical Algorithms, 2nd Edition. Willmann-Bell, Inc., Richmond, VA.
Seidelmann, P.K. 2006. Explanatory Supplement to the Astronomical Almanac. University
Science Books, Sausalito, CA.

Seidelmann, P.K. 2006. Explanatory Supplement to the Astronomical Almanac. University Science Books, Sausalito, CA. 
This function is an R implementation of: 

Jacobson L, Seaver A, Tang J. 2011. AstroCalc4R: software to calculate solar zenith angle; 
time at sunrise, local noon and sunset; and photosynthetically available radiation based on date, 
time and location. US Dept Commer, Northeast Fish Sci Cent Ref Doc. 11-14; 10 p.  
Available from: National Marine Fisheries Service, 166 Water Street, Woods Hole, MA 02543-1026, 
or online at \url{http://www.nefsc.noaa.gov/nefsc/publications/}.
Note that the documentation may be easiest to obtain at
\url{http://www.nefsc.noaa.gov/publications/crd/crd1114/index.html} 
}
\author{
Larry Jacobson, Alan Seaver, and Jiashen Tang
NOAA National Marine Fisheries Service
Northeast Fisheries Science Center, 166 Water St., Woods Hole, MA 02543

\email{Larry.Jacobson@noaa.gov}
}
\examples{
astrocalc4r(day=12,month=9,year=2000,hour=12,timezone=-5,lat=40.9,lon=-110)
}
\keyword{MISC}
