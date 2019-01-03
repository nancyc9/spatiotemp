
# This code is to clean the timeseries of data for all sites.
# NOTE: This function does not compute time trends, that happens later
# This function ignores issues related to LOD and analyst flags.
# The mesa data is assumed to have a column called "site_type" which is either
# F for fixed, C for comco, or O for home.
# Both aqsData and mesaData are assumed to have a column called location_id.
# aqsConcColumn and mesaConcColumn are the names of the concentration columns
# in the two datasets.
# Both aqsData and mesaData are assumed to have a column called intended_wednesday,
# which is the date that will be used for each observation record, and a
# coloumn called location_id, which identifies the measurement location.
# npls is the number of pls components that will be used.
timeseries <- function(aqsData, mesaData, aqsConcColumn, mesaConcColumn, npls) {


  ########################### Checking NAs

  if (!(aqsConcColumn %in% colnames(aqsData))) {
    stop("the concentration column is missing from the aqs data")
  }
  if (!(mesaConcColumn %in% colnames(mesaData))) {
    stop("the concentration column is missing from the mesa data")
  }

  # Check to make sure all dates are wednesdays.
  if (any(weekdays(aqsData$intended_wednesday) != "Sunday")) {
    stop("aqsData$intended_wednesday must all be Sundays")
  }
  if (any(weekdays(mesaData$intended_wednesday) != "Sunday")) {
    print(weekdays(mesaData$intended_wednesday[which(weekdays(mesaData$intended_wednesday) != "Sunday")]))
    stop("mesaData$intended_wednesday must all be Sundays")
  }

  # harmonize names
  aqsData$concentration <- aqsData[,aqsConcColumn]
  mesaData$concentration <- mesaData[,mesaConcColumn]
  aqsData$ID <- aqsData$native_id
  mesaData$ID <- mesaData$native_id
  mesaData$type <- "unknown"
  mesaData[which(mesaData$site_type == "F"), ]$type <- "FIXED"
  if (any(mesaData$site_type == "C")) {
    mesaData[which(mesaData$site_type == "C"), ]$type <- "COMCO"
  }
  if (any(mesaData$site_type == "O")) {
    mesaData[which(mesaData$site_type == "O"), ]$type <- "HOME"
  }
  print(paste("timeseries: failed to label ", dim(mesaData[which(mesaData$type=="unknown"),])[1], " mesa sites.", sep=""))
  aqsData$type <- "AQS"

  # Check for NAs in concentration values
  if (sum(is.na(aqsData$concentration)) != 0) {
    warning(paste("removing ",sum(is.na(aqsData$concentration))," NA values from aqs concentration data", sep=""))
    aqsData <- aqsData[!is.na(aqsData$concentration), ]
  }
  if (sum(is.na(mesaData$concentration)) != 0) {
    warning(paste("removing ",sum(is.na(mesaData$concentration))," NA values from mesa concentration data", sep=""))
    mesaData <- mesaData[!is.na(mesaData$concentration), ]
  }

  # Check for zeros or negative values in concentration values
  if (sum(aqsData$concentration<=0) != 0) {
    warning("removed a <=0 value in aqs concentration data")
    aqsData <- aqsData[aqsData$concentration>0, ]
  }
  if (sum(mesaData$concentration<=0) != 0) {
    warning("removed a <=0 value in mesa concentration data")
    mesaData <- mesaData[mesaData$concentration>0, ]
  }

  # remove any non-comco and non-home measurement locations where the number of observations
  # is <= npls, to avoid problems with creating pls time trend regressions.
  # COMCO and HOME data is not used in pls time trend regressions.
  x <- data.frame(table(mesaData$ID))
  y <- x[which(x$Freq > npls),]
  mesaData <- mesaData[which((mesaData$native_id %in% y$Var1) | mesaData$type == "COMCO" | mesaData$type == "HOME"), ]
  x2 <- data.frame(table(mesaData$ID))
  if (dim(x)[1] != dim(x2)[1]) {
   warning(paste("removed ", dim(x)[1] - dim(x2)[1], " mesa measurement locations ",
   "where the number of measurements was <= npls (", npls, ")", sep=""))
  }

  x <- data.frame(table(aqsData$ID))
  y <- x[which(x$Freq > npls),]
  aqsData <- aqsData[which((aqsData$native_id %in% y$Var1) | aqsData$type == "COMCO" | aqsData$type == "HOME"), ]
  x2 <- data.frame(table(aqsData$ID))
  if (dim(x)[1] != dim(x2)[1]) {
   warning(paste("removed ", dim(x)[1] - dim(x2)[1], " aqs measurement locations ",
   "where the number of measurements was <= npls (", npls, ")", sep=""))
  }

  # combine AQS 2-week averaged data with mesa 2-week average data #
  allData <- rbind(
    mesaData[, c("ID", "intended_wednesday", "concentration", "type")],
    aqsData[, c("ID", "intended_wednesday", "concentration", "type")],
    make.row.names = FALSE)



  time_series <- createDataMatrix(obs=allData$concentration, date=allData$intended_wednesday, ID=allData$ID)

  out <- list(data = allData, matrix = time_series)
  return(out)
}
