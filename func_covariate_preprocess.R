
#==============================================================================================================================================================================
# 'covariate.preprocess'
#
# Input arguments:
# * covars.mon: region-specific covariate data at monitoirng sites
# * covars.sub: region-specific covariate data at cohort addresses
#     covariate data for moniotirng and cohort locations should
#        # be region-specific
#        # be selected for the defined areas
#        # be separate by monitoring and cohort locations
#        # include the same covariates between montioring and cohort data
#        # have no missing values for any variables (residual oil for NY, street canyon for NY and IL, and CALINE alternative for CA)
#        # have the column, 'type', for type of site including "AQS", "FIXED", "HOME", "COMCO" or "P"
# * region: study region for monitoirng and cohort data ('CA', 'NY', and 'IL' should be used when analysts want to include CALINE alternative, residual oil, or street canyon)
# * exclude.cat.canyon: preprocessing 2-0. whether or not to exclude non-numeric street canyon variable (canyon_road_type) (default is TRUE)
# * exclude.cityhall:   preprocessing 2-0. whether or not to exclude distances to main and local city halls (default is TRUE)
# * exclude.census:     preprocessing 2-0. whether or not to exclude census variables (default is TRUE)
# * exclude.old.LUR:    preprocessing 2-0. whether or not to exclude USGS land use variables constructed in 1970-80 (default is TRUE)
# * exclude.relative.elevation: preprocessing 2-0. whether or not to exclude relative elevation (default is FALSE)
#
# Output objects:
# * covars.mon: preprocessed region-specific covariate data at monitoring astes
# * covars.sub: preprocessed region-specific covarirate data at cohort addresses
# * exclude.include.vars: initial, excluded, and included variables
#       * common.value: preprocessing 2-1. list of excluded variables which have less than 20% of being different from the most common value in monitoring data
#       * low.landuse:  preprocessing 2-2. list of excluded varlabies with maximum land use variables less than 10 in monitoring data
#       * sd.ratio:     preprocessing 2-3. list of excluded variables with SD of cohort data greater than 5 times of SD of monitoring data
#       * outlier:      preprocessing 2-4. list of excluded variables with outliers more than 2% in monitoring and cohort data combined
#       * initial:      list of initial variables after preprocessing 2-0
#       * include:      list of final variables after preprocessing 2
# * vars.count: numbers of excluded variables for each preprocessing step
#==============================================================================================================================================================================


covariate.preprocess <- function(covars.mon, covars.sub, region,
    exclude.cat.canyon = TRUE, exclude.cityhall = TRUE, exclude.census = TRUE,
    exclude.old.LUR = TRUE, exclude.relative.elevation = FALSE )
{

  if (dim(covars.mon)[1]==0) {
     stop("There are zero rows of monitor data")
  }
  if (dim(covars.sub)[1]==0) {
     stop("There are zero rows of subject data")
  }

# 0. Data prepration
    if(region!="NY" & region!="IL"){
        oil <- grepl("oil", names(covars.mon))
        canyon <- grepl("canyon", names(covars.mon))                    #canyon: street canyon variables for NY and Chicago
        bus <- grepl("bus", names(covars.mon))
        covars.mon <- covars.mon[,!oil & !canyon & !bus]
        covars.sub <- covars.sub[,!oil & !canyon & !bus]
    }
    if(region=="IL"){
        oil <- grepl("oil", names(covars.mon))                          #oil: residual oil variables for NY
        bus <- grepl("bus", names(covars.mon))
        covars.mon <- covars.mon[,!oil & !bus]
        covars.sub <- covars.sub[,!oil & !bus]
    }
    if(region!="CA"){
        caline.alt <- grepl('^calinemod_alt_lt',names(covars.mon))      #caline.alt: caline alternative variables for LA
        covars.mon <- covars.mon[,!caline.alt]
        covars.sub <- covars.sub[,!caline.alt]
    }

    covars.sub <- covars.sub[,names(covars.mon)]
    covars.all <- rbind(covars.mon, covars.sub)

# 1. Preprocessing 1: recoding
# 1.1. compute min distance to any roads (a1, a2 and m3)
    covars.a123D <- combine_a123_m(covars.all, removeOrig=FALSE)
# 1.2. comute and small roads (a2 and a3) and remove original variables
    covars.a23D  <- combine_a23_m(covars.a123D, removeOrig=TRUE)
# 1.3. comute and small roads (a2 and a3) and remove original variables
    covars.a23L  <- combine_a23_ll(covars.a23D, removeOrig=TRUE)

# 1.4. natural log transform CALINE, emission, and distance to sources (roads, any & large airports, coast, SML ports, commercial area, railrod/yard, oil, and cityhall)
#      truncate distance variables at 10 m
#      remove original variables
    covars.cal      <- log_transform_caline(covars.a23L, removeOrig=TRUE)
    covars.em       <- log_transform_emission(covars.cal, removeOrig=TRUE)
    covars.process1 <- log_transform_distance(covars.em, lowerBound=10, removeOrig=TRUE)

# 2. Preprocessing 2: variable exclusion
    covars.process1.mon <- covars.process1[covars.process1$type!="P",]
    covars.process1.sub <- covars.process1[covars.process1$type=="P",]

# 2.0. exclude distances to road types for street canyons, city halls, old LUR, census, or relative elevation variables
    gis.vars <- c('m_to|pop|ll|ndvi|imp|em|calinemod|rlu|lu|elev|tl|intersect|oil|canyon')
    initial.vars <- names(covars.process1)[grepl(gis.vars, names(covars.process1))]

    cityhall <- names(covars.process1)[grep('cityhall',names(covars.process1))]
    old.LUR <- names(covars.process1)[grep('^lu',names(covars.process1))]
    census <- names(covars.process1)[grep('^tr|^bg|^bk',names(covars.process1))]
    relative.elevation <- names(covars.process1)[grep('elev_1k|elev_5k',names(covars.process1))]

    if(exclude.cat.canyon) {
        initial.vars <- initial.vars[!(initial.vars %in% 'canyon_road_type')]
    }
    if(exclude.cityhall) {
        initial.vars <- initial.vars[!(initial.vars %in% cityhall)]
    }
    if(exclude.old.LUR) {
        initial.vars <- initial.vars[!(initial.vars %in% old.LUR)]
    }
    if(exclude.census) {
        initial.vars <- initial.vars[!(initial.vars %in% census)]
    }
    if(exclude.relative.elevation) {
        initial.vars <- initial.vars[!(initial.vars %in% relative.elevation)]
    }

# 2.1 exclude variables with values less than 20% of being different from most common values
    # For this step we only use fixed and AQS sites, because when we normalize the covariates
    # later we do it based on the standard devaiation of the fixed and AQS sites.
    fixed_and_aqs_indices <- which(covars.process1.mon$type == "AQS" | covars.process1.mon$type == "FIXED")

    exclude.common.value  <- fail_most_common_value(covars.process1.mon[fixed_and_aqs_indices, ], initial.vars, thres=0.2)

# 2.2 exclude variables max landuse variable < 10%
    exclude.low.landuse   <- fail_low_landuse(covars.process1.mon, initial.vars, lowValue=10)

# 2.3 exclude variables with sd of cohort > 5 times sd of monitoring data
    exclude.sd.ratio      <- fail_sd_ratio(covars.process1.mon, covars.process1.sub, initial.vars, thres=5)

# 2.4 exclude variables more than 2% outliers in monitor sites and ppt sites combined
    exclude.outlier       <- fail_outlier_check(covars.process1.mon, covars.process1.sub, initial.vars, thres=0.02)
    exclude.vars <- list(exclude.common.value, exclude.low.landuse, exclude.sd.ratio, exclude.outlier)
    names(exclude.vars) <- c('common.value','low.landuse','sd.ratio','outlier')
    include.vars <- initial.vars[!(initial.vars %in% unique(unlist(exclude.vars)))]

    covars.process2.mon <- covars.process1.mon[,include.vars]
    covars.process2.sub <- covars.process1.sub[,include.vars]

# 3. Output extraction
    exclude.include.vars <- list(exclude.common.value, exclude.low.landuse,
        exclude.sd.ratio, exclude.outlier, initial.vars, include.vars)
    names(exclude.include.vars) <- c(names(exclude.vars), 'initial', 'include')
    vars.count <- c(length(initial.vars), sapply(exclude.vars,length), length(include.vars))
    names(vars.count) <- c('Initial',names(exclude.vars),'Include')

    covars <- list(covars.mon=covars.process2.mon, covars.sub=covars.process2.sub,
        exclude.include.vars=exclude.include.vars, vars.count=vars.count)
    return(covars)
}




#======================================================================================================================
# Sub-functions
# 1-1-1) combine_a123_m:         min distance to any roads (a1,a2,a3)
# 1-1-2) combine_a23_m:          min of a2,a3 distance
# 1-1-3) combine_a23_ll:         sum of a2 and a3 roads in buffers
# 1-2-1) log_transform_caline:   natural log transform CALINE variables
# 1-2-2) log_transform_emission: natural log transform emission variables
# 1-2-3) log_transform_distance: natural log transform distance variables and truncate at 10 meters
# 2-1) fail_most_common.value:   variables with less than <20> % being differing from the most common value
# 2-2) fail_low_landuse:         land use variables (old and new) whose max value is less than <10> % in monitoring data
# 2-3) fail_sd_ratio:            variables whose SD for cohort data is <5> times greater than SD for monitoirng data
# 2-4) fail_outlier_check:       variables with more than <2> % outliers (Z score >5)
#=======================================================================================================================

combine_a123_m <- function(all.data, removeOrig=FALSE)
{
  ma1a2a3.index <- grepl("^m_to_a[123]$", colnames(all.data))
  newcol.index <- 1 + ncol(all.data)
  all.data[, newcol.index] <- apply(all.data[, ma1a2a3.index], 1, min)
  colnames(all.data)[newcol.index] <- 'm_to_a123'
  if (removeOrig) all.data <- all.data[, !ma1a2a3.index]
  return(all.data)
}


combine_a23_m <- function(all.data, removeOrig=FALSE)
{
  ma2a3.index <- grepl("^m_to_a[23]$", colnames(all.data))
  newcol.index <- 1 + ncol(all.data)
  all.data[, newcol.index] <- apply(all.data[, ma2a3.index], 1, min)
  colnames(all.data)[newcol.index] <- 'm_to_a23'
  if (removeOrig) all.data <- all.data[, !ma2a3.index]
  return(all.data)
}


combine_a23_ll <- function(all.data, removeOrig=FALSE)
{
  a2.vars <- grep("ll_a2", colnames(all.data))
  if (length(a2.vars) == 0) { return(all.data)}
  for (i in a2.vars)
  {
    newcol.index <- 1 + ncol(all.data)
    a3.var <- grep(gsub('a2','a3',colnames(all.data)[i]), colnames(all.data))
    all.data[, newcol.index] <- all.data[, i] + all.data[, a3.var]
    colnames(all.data)[newcol.index] <- paste("ll_a23_", strsplit(colnames(all.data)[a3.var], '_')[[1]][3], sep="")
  }
  ll.vars <- grep("ll_a[^1]_s", colnames(all.data))
  if (removeOrig) all.data <- all.data[, -ll.vars]
  return(all.data)
}


log_transform_caline <- function(all.data, removeOrig=FALSE)
{
  caline.vars <- grep("caline", colnames(all.data))
  if (length(caline.vars) == 0) { return(all.data)}
  new.varnames <- c()
  for (i in caline.vars)
  {
    newcol.index <- 1 + ncol(all.data)
    all.data[, newcol.index] <- log(all.data[, i] + 0.1)
    colnames(all.data)[newcol.index] <- paste('log_', colnames(all.data)[i], sep='')
  }
  if (length(caline.vars)>0 & removeOrig) all.data <- all.data[, -caline.vars]
  return (all.data)
}


log_transform_emission <- function(all.data, removeOrig=FALSE)
{
  em.vars <- grep("^em_", colnames(all.data))
  if (length(em.vars) == 0) { return(all.data)}
  new.varnames <- c()
  for (i in em.vars)
  {
    newcol.index <- 1 + length(colnames(all.data))
    all.data[, newcol.index] <- log(all.data[, i] + 0.1)
    colnames(all.data)[newcol.index] <- paste('log_', colnames(all.data)[i], sep='')
  }
  if (removeOrig) all.data <- all.data[, -em.vars]
  return (all.data)
}


log_transform_distance <- function(all.data, lowerBound=10, removeOrig=FALSE)
{
  distance.vars <- grep("^m_to_", colnames(all.data))
  if (length(distance.vars) == 0) { return(all.data)}
  new.varnames <- c()
  for (i in distance.vars)
  {
    newcol.index <- 1 + ncol(all.data)
    all.data[, newcol.index] <- log( sapply( all.data[, i], function(x) { max(lowerBound, x) } ) )
    colnames(all.data)[newcol.index] <- paste('log_', colnames(all.data)[i], sep='')
  }
  if (removeOrig) all.data <- all.data[, -distance.vars]
  return (all.data)
}


fail_most_common_value <- function(mon.data, vars.all, thres=0.2)
{
  thres <- dim(mon.data[,vars.all])[1]*thres
  fail <- apply( mon.data[,vars.all], 2, function(x) {
    tmp <- split(x,x)
    most.common.value <- tmp[[which.max(sapply(tmp, length))]][1]
    sum(x != most.common.value, na.rm=T) } ) < thres
  fail <- names(mon.data[,vars.all])[fail]
  return(fail)
}


fail_low_landuse <- function(mon.data, vars.all, lowValue=10)
{
  lu.vars <- grep("^rlu_|^lu", vars.all, value=T)
  fail <- sapply(lu.vars, function(x) return (max(mon.data[, x]) < lowValue))
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}


fail_sd_ratio <- function(mon.data, cohort.data, vars.all, thres=5)
{
  fail <- c()
  for (i in vars.all)  {
    mon.sd <- sd(mon.data[, i], na.rm=TRUE)
    cohort.sd <- sd(cohort.data[, i], na.rm=TRUE)
    if (cohort.sd > thres * mon.sd | all(is.na(mon.sd)))
      fail <- c(fail, i)
  }
  return (fail)
}


fail_outlier_check <- function(mon.data, cohort.data, vars.all, thres=0.02)
{
  all.data = rbind(mon.data[,vars.all], cohort.data[,vars.all])

  fail <- sapply( vars.all, function(x)
    return (  sum(abs(scale(all.data[, x]))> 5) > nrow(all.data)*thres  ) )
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}
