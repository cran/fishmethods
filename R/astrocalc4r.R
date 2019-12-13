astrocalc4r=function (day, month, year, hour, timezone, lat, lon, withinput = FALSE, 
          seaorland = "maritime", acknowledgment = FALSE) 
{
  if (acknowledgment) {
    cat("\n", "---------------------------------------------------------")
    cat("\n", "                AstroCalcPureR Version 2.3")
    cat("\n", "Documentation: Jacobson L, Seaver A, Tang J. 2011. AstroCalc4R:")
    cat("\n", "software to calculate solar zenith angle; time at sunrise,")
    cat("\n", "local noon and sunset; and photosynthetically available")
    cat("\n", "radiation based on date, time and location. US Dept Commer,")
    cat("\n", "Northeast Fish Sci Cent Ref Doc. 11-14; 10 p. Available from:")
    cat("\n", "National Marine Fisheries Service, 166 Water Street, ")
    cat("\n", "Woods Hole, MA 02543-1026, or online at")
    cat("\n", "http://nefsc.noaa.gov/publications/")
    cat("\n \n", "Available in fishmethods library.  Contact the fishmethods")
    cat("\n" , "administrator or Larry Jacobson (NOAA, National Marine")
    cat("\n", "Fisheries Service-retired) at larryjacobson6@gmail.com")
    cat("\n", "for assitance.")
    cat("\n\n", "Useage:")
    cat("\n", "    AstroCalcPureR(day,month,year,hour,timezone,")
    cat("\n", "                   lat,lon,withinput=F,")
    cat("\n", "                   seaorland='maritime',")
    cat("\n", "                   acknowledgment=TRUE)")
    cat("\n\n", "HINT: set acknowledgment=FALSE to avoid this message")
    cat("\n", "---------------------------------------------------------", 
        "\n")
  }
  ###
  # Modification history:
  # Added traps to ensure that timezone and longitude have the same sign if
  # time is not UTC (timezone==0).  In particular, if time is not UTC
  # then both should be negative in western hemisphere or both positive in 
  # the eastern hemisphere.  This is a bit tricky because
  # the boundaries for time zones are not based strictly on latitude.  The 
  # solution here is to require UTC (universal stardard time inputs) input data 
  # (argument timezone=0) if  there is a conflict.  In other words, require
  # the user to convert input time to UTC and use timezone=0 if necessary.
  #
  # Some minor changes to error messages.
  #
  # Version number in optional acknoledgement changed from to 2.2 to 2.3.
  #                        - Larry Jacobson, January 13, 2019
  #
  # timezone=0 if a sign conflict
  options(digits = 9)
  deg2rad <- pi/180
  null.c <- function(x) return(sum(is.null(x)))
  if (sum(null.c(day), null.c(month), null.c(year), null.c(hour), 
          null.c(timezone), null.c(lat), null.c(lon)) > 0) 
    stop("\n Null or missing required data vector for day, month, year, timezone, lat or lon \n")
  if ((length(day) != length(month)) | (length(month) != length(year)) | 
      (length(year) != length(hour)) | (length(hour) != length(timezone)) | 
      (length(timezone) != length(lat)) | (length(lat) != 
                                           length(lon))) 
    stop("\n Input vectors are not the same length \n")
  times <- length(day)
  na.c <- function(x) return(sum(is.na(x)))
  if (sum(na.c(day), na.c(month), na.c(year), na.c(hour), 
          na.c(timezone), na.c(lat), na.c(lon)) > 0) 
    stop("\n NA values in input data \n")
  logic1 <- year < 0
  if (sum(logic1) > 0) 
    stop(cat("\n Error in year at rows:", (1:times)[logic1], 
             " \n\n"))
  is.leap <- function(x) return((((x%%4 == 0) & (x%%100 != 
                                                   0))) | (x%%400 == 0))
  date.list <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 
                 31)
  logic1 <- abs(month - 6) > 6
  if (sum(logic1) > 0) 
    stop(cat("\n Error in month at rows:", (1:times)[logic1], 
             " \n\n"))
  logic1 <- day > (date.list[month] + is.leap(year) * (month == 
                                                         2))
  logic2 <- day <= 0
  if ((sum(logic1) > 0) | (sum(logic2) > 0)) 
    stop(cat("\n Incorrect month-day-year combination at rows: ", 
             (1:times)[logic1 | logic2], " \n\n"))
  logic1 <- abs(hour - 12) > 12
  if (sum(logic1) > 0) 
    stop(cat("\n Error in hour at rows:", (1:times)[logic1], 
             " \n\n"))
  logic1 <- abs(timezone) > 12
  if (sum(logic1) > 0) 
    stop(cat("\n Error in time zone at rows:", 
             (1:times)[logic1], " \n\n"))
  logic1 <- abs(lon) > 180
  if (sum(logic1) > 0) 
    stop(cat("\n Error in longitude at rows:", 
             (1:times)[logic1], " \n\n"))
  logic1 <- abs(lat) > 90
  if (sum(logic1) > 0) 
    stop(cat("\n Error in latitude at rows:", 
             (1:times)[logic1], " \n\n"))
  logic1 <- sign(lon) == sign(timezone)
  logic2 <- timezone == 0
  logic3 <- !(logic1 | logic2)
#  print(c(logic1=logic1,logic2=logic2,logic3))
  if(sum(logic3) !=0) stop(cat("\n \n Arguments longitude and timezone must have the same sign if input time is",
                                    "\n not UTC (timezone != 0).  In particular, if timezone !=0, both lon and timezone must",
                                    "\n be negative for locations in western hemisphere and positive for locations in the",
                                    "\n eastern hemisphere.  Check and fix input data if warranted. If data are correct",
                                    "\n then convert input time (argument hour) to UTC and use timezone=zero.", 
                                    "\n This problem  occurs ",sum(logic3)," times at rows: ",(1:times)[logic3],"\n\n"))
  
  JulianDay <- function(xday, xmonth, xyear) {
    mm <- xmonth
    xmonth[mm <= 2] <- xmonth[mm <= 2] + 12
    xyear[mm <= 2] <- xyear[mm <= 2] - 1
    xa <- floor(xyear/100)
    xb <- 2 - xa + floor(xa/4)
    jd <- floor(365.25 * (xyear + 4716)) + floor(30.6001 * 
                                                   (xmonth + 1)) + xday + xb - 1524.5
    return(jd)
  }
  daymonth <- function(mth, yr) {
    day[is.leap(yr)] <- c(31, 29, 31, 30, 31, 30, 31, 31, 
                          30, 31, 30, 31)[mth[is.leap(yr)]]
    day[!is.leap(yr)] <- c(31, 28, 31, 30, 31, 30, 31, 31, 
                           30, 31, 30, 31)[mth[!is.leap(yr)]]
    return(day)
  }
  parcalc <- function(zenith, setting = seaorland) {
    I0 <- 531.2
    V <- 23
    uv <- 1.4
    u0 <- 0.34
    r <- 0.05
    d <- 1
    if (!setting %in% c("maritime", "continental")) 
      stop("setting value is neither 'maritime' nor 'continental'!")
    if (setting == "maritime") {
      a <- 0.068
      b <- 0.379
      a1 <- 0.117
      b1 <- 0.493
      av <- 0.002
      bv <- 0.87
      a0 <- 0.052
      b0 <- 0.99
    }
    else if (setting == "continental") {
      a <- 0.078
      b <- 0.882
      a1 <- 0.123
      b1 <- 0.594
      av <- 0.002
      bv <- 0.87
      a0 <- 0.052
      b0 <- 0.99
    }
    zrad <- zenith * deg2rad
    x1 <- uv/cos(zrad)
    xx <- exp(-av * x1^bv)
    x2 <- u0/cos(zrad)
    xxx <- exp(-a0 * x2^b0)
    xa <- a + b/V
    xb <- d - r * (a1 + b1/V)
    par <- I0 * cos(zrad) * exp(-xa/cos(zrad))/xb * xx * 
      xxx
    par[zenith > 89.9999] <- 0
    return(par)
  }
  output <- as.data.frame(matrix(nrow = 0, ncol = 9))
  names(output) <- c("noon", "sunrise", "sunset", "azimuth", 
                     "zenith", "eqtime", "declin", "daylight", "PAR")
  hourtemp <- hour - timezone
  hour <- ifelse(hourtemp > 24, hourtemp - 24, hourtemp)
  change_day <- !(hour == hourtemp)
  dm <- daymonth(month, year)
  daytemp <- day
  daytemp[change_day] <- ifelse((day[change_day] < dm[change_day]), 
                                day[change_day] + 1, 1)
  change_month <- abs(day - daytemp) > 1
  monthtemp <- month
  monthtemp[change_month] <- ifelse(month[change_month] < 
                                      12, month[change_month] + 1, 1)
  change_year <- abs(month - monthtemp) > 1
  yeartemp <- year
  yeartemp[change_year] <- year[change_year] + 1
  xy <- yeartemp
  xm <- monthtemp
  xd <- daytemp + hourtemp/24
  jd <- JulianDay(xd, xm, xy) * 100/100
  jc <- (jd - 2451545)/36525
  xx <- 280.46646 + 36000.76983 * jc + 0.0003032 * jc^2
  gmls <- xx%%360
  xx <- 357.52911 + 35999.05029 * jc - 0.0001537 * jc^2
  gmas <- xx%%360
  eeo <- 0.016708634 - 4.2037e-05 * jc - 1.267e-07 * jc^2
  scx <- (1.914602 - 0.004817 * jc - 1.4e-05 * jc^2) * sin(gmas * 
                                                             deg2rad) + (0.019993 - 0.000101 * jc) * sin(2 * gmas * 
                                                                                                           deg2rad) + 0.000289 * sin(3 * gmas * deg2rad)
  Stl <- gmls + scx
  Sta <- gmas + scx
  srv <- 1.000001018 * (1 - eeo^2)/(1 + eeo * cos(Sta * deg2rad))
  omega <- 125.04 - 1934.136 * jc
  lambda <- Stl - 0.00569 - 0.00478 * sin(omega * deg2rad)
  epsilon <- (23 + 26/60 + 21.448/60^2) - (46.815/60^2) * 
    jc - (0.00059/60^2) * jc^2 + (0.001813/60^2) * jc^3
  oblx <- 0.00256 * cos(omega * deg2rad)
  epsilon <- epsilon + oblx
  alpha <- atan2(cos(epsilon * deg2rad) * sin(lambda * deg2rad), 
                 cos(lambda * deg2rad))/deg2rad
  declin <- asin(sin(epsilon * deg2rad) * sin(lambda * deg2rad))/deg2rad
  y <- tan(epsilon * deg2rad/2)^2
  eqtime <- (y * sin(2 * gmls * deg2rad) - 2 * eeo * sin(gmas * 
                                                           deg2rad) + 4 * eeo * y * sin(gmas * deg2rad) * cos(2 * 
                                                                                                                gmls * deg2rad) - y^2 * sin(4 * gmls * deg2rad)/2 - 
               5/4 * eeo^2 * sin(2 * gmas * deg2rad))/deg2rad * 4
  h0 <- -0.8333 * deg2rad
  phi <- lat * deg2rad
  hangle <- acos((sin(h0) - sin(declin * deg2rad) * sin(phi))/cos(declin * 
                                                                    deg2rad)/cos(phi))/deg2rad
  noon <- (720 - 4 * lon + timezone * 60 - eqtime)/1440
  sunrise <- (noon * 1440 - hangle * 4)/1440 * 24
  sunset <- (noon * 1440 + hangle * 4)/1440 * 24
  noon <- noon * 24
  daylight <- hangle * 8
  tst <- (hourtemp * 60 + eqtime + 4 * lon)%%1440
  tsa <- ifelse(tst < 0, tst/4 + 180, tst/4 - 180)
  zenith <- 90 - asin(sin(lat * deg2rad) * sin(declin * deg2rad) + 
                        cos(lat * deg2rad) * cos(declin * deg2rad) * cos(tsa * 
                                                                           deg2rad))/deg2rad
  azimuth <- acos((sin(lat * deg2rad) * sin((90 - zenith) * 
                                              deg2rad) - sin(declin * deg2rad))/cos(lat * deg2rad)/cos((90 - 
                                                                                                          zenith) * deg2rad))/deg2rad + 180
  azimuth <- ifelse(tsa > 0, azimuth%%360, 360 - azimuth%%360)
  daylight <- daylight/60
  PAR <- parcalc(zenith)
  if (any(is.nan(sunrise))) {
    message(paste("Warning: Polar day/night (daylength 0 or 24 hrs) at record(s):", 
                  (1:times)[is.nan(sunrise)], "\n Check input data (i.e. latitude)?"))
    daylight <- ifelse(PAR > 0, 24, 0)
  }
  output <- rbind(output, data.frame(noon = noon, sunrise = sunrise, 
                                     sunset = sunset, azimuth = azimuth, zenith = zenith, 
                                     eqtime = eqtime, declin = declin, daylight = daylight, 
                                     PAR = PAR))
  if (withinput) 
    return(cbind(data.frame(tzone = timezone, day = day, 
                            month = month, year = year, hhour = hour, xlat = lat, 
                            xlon = lon), output))
  else return(output)
}
