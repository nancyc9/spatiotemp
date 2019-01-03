
# This function prepares the covariate data.
# covariates is the input data frame. It is expected to have columns named
# region, location_id, longitude, latitude
# obs is a data list of the type output by the "timeseries" function.
# regionName is the name of the region being processed. It will be matched to the
# "region" column in the covariates data frame.
# This function makes up random participant data to satisfy the requirements of
# other functions. This function will need to be edited in order to be used
# with actual participant data.
covariatePrepare <- function(covariates, obs, regionName) {
  print(paste("covariatePrepare: input observation data has ",sum(!is.na(obs$matrix))," measurements", sep=""))

  covariates$type <- "unknown"

  aqs_idx <- which(obs$data$type == "AQS")
  aqs_indices <- covariates$native_id %in% obs$data[aqs_idx, "ID"]
  covariates[aqs_indices, ]$type <- "AQS"

  fixed_obs_idx <- which(obs$data$type == "FIXED")
  fixed_indices <- covariates$native_id %in% obs$data[fixed_obs_idx, "ID"]
  if (sum(fixed_indices) == 0) {
    print(obs$data[fixed_obs_idx, "ID"])
    print(covariates$native_id)
    stop("There do not exist any covariates with IDs matching the fixed measurements")
  }
  covariates[fixed_indices, ]$type <- "FIXED"

  if (any(obs$data$type == "COMCO")) {
    comco_obs_idx <- which(obs$data$type == "COMCO")
    comco_indices <- covariates$native_id %in% obs$data[comco_obs_idx, "ID"]
    covariates[comco_indices, ]$type <- "COMCO"
  }
  if (any(obs$data$type == "HOME")) {
    home_obs_idx <- which(obs$data$type == "HOME")
    home_indices <- covariates$native_id %in% obs$data[home_obs_idx, "ID"]
    covariates[home_indices, ]$type <- "HOME"
  }

  # output data
  monitor_indices <- which(covariates$type == "AQS" | covariates$type ==
      "FIXED" | covariates$type == "COMCO" | covariates$type == "HOME")
  covars.monitors <- covariates[monitor_indices, ]
  print(paste("covariatePrepare: Keeping covariates for ", dim(covars.monitors)[1], " of ", dim(covariates)[1], " locations.", sep=""))
  if (dim(covars.monitors)[1] < dim(covariates)[1]) {
    print("dropping covariate sites that are not AQS, FIXED, COMCO, or HOME: ")
    print(covariates[which(covariates$type != "AQS" & covariates$type !=
        "FIXED" & covariates$type != "COMCO" & covariates$type != "HOME"),"native_id"])
  }
  rownames(covars.monitors) <- as.character(covars.monitors$native_id)
  # Split out 25% of the data to consider as participants.  This is necessary
  # because the functions require participant data but we don't have any.
  smp_size <- floor(0.75 * nrow(covars.monitors))
  train_ind <- sample(seq_len(nrow(covars.monitors)), size = smp_size)
  covars.ppt <- covars.monitors[-train_ind, ]
  covars.ppt$type <- "P"

  d <- list() # Data holder for this region
  d$map <- regionName
  # Figure out which rows are for the current region.
  region_indices <- which(covars.monitors$region == regionName)
  print(paste("covariatePrepare: keeping the ",length(region_indices), " of ", dim(covars.monitors)[1],
    " total covariate locations that are in the ", regionName, " region.", sep=""))

  d$covars.monitors <- covars.monitors[region_indices, ]
  # use the same covars.ppt for each region because it doesn't matter
  d$covars.ppt <- covars.ppt


  # Creating the st_matrix
  region_monitors <- colnames(obs$matrix) %in% d$covars.monitors$native_id
  print(paste("Keeping observation data for ",sum(region_monitors, na.rm=TRUE), " of ",dim(obs$matrix)[2], " monitor sites.", sep=""))
  if (sum(region_monitors) < dim(obs$matrix)[2]) {
    print("missing covariate data: dropping measurements for ids: ")
    print(colnames(obs$matrix)[which(!region_monitors)])
  }
  d$st_matrix <- obs$matrix[, region_monitors]

  # Reorder so the locations match the covariates
  d$st_matrix <- d$st_matrix[,match(d$covars.monitors$native_id, colnames(d$st_matrix))]

  # longitude and latitude
  d$lon_lat <- d$covars.monitors[, c("longitude", "latitude")]
  d$lon_lat_ppt <- d$covars.ppt[, c("longitude", "latitude")]
  # lambert coords for modeling monitors and for ppt locations too
  d$lambert.monitors <- d$covars.monitors[, c("lambert_x", "lambert_y")]
  d$lambert.ppt <- d$covars.ppt[, c("lambert_x", "lambert_y")]
  d$type.monitors <- d$covars.monitors$type
  d$type.ppt <- d$covars.ppt$type
  d$native_id.monitors <- d$covars.monitors$native_id
  d$native_id.ppt <- d$covars.ppt$native_id

  # Create the 'covars' object.
  # Includes preprocessing of covariates and exclusion criteria
  covars.temp <- covariate.preprocess(covars.mon = d$covars.monitors,
      covars.sub = d$covars.ppt, region = d$map)

  d$covars.monitors <- covars.temp$covars.mon
  d$covars.ppt <- covars.temp$covars.sub

  d$vars.count <- covars.temp$vars.count
  d$exclude.include.vars <- covars.temp$exclude.include.vars

  # Drop covariates that have all the same value in the AQS and Fixed sites.
  # This is found by looking for covariates with sd = 0. They pose problems
  # for the scaling done in PLS.
  d$aqs_and_fixed_indices <- which(d$type.monitors == "AQS" | d$type.monitors == "FIXED")
  bad_covar <- which(apply(d$covars.monitors[d$aqs_and_fixed_indices, ],
     2, sd, na.rm = T) < 0.001)
  if (length(bad_covar) > 0) {
      d$covars.monitors <- d$covars.monitors[, -bad_covar]
      d$covars.ppt <- d$covars.ppt[, -bad_covar]
  }
  d$exclude.include.vars$all_same_in_fixed <- names(bad_covar)
  print(paste("covariatePrepare: output observation data has ",sum(!is.na(d$st_matrix))," measurements", sep=""))
  return(d)
}
