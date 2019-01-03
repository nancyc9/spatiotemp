# model_fit_params creates, saves, and returns and ST model estimate, where
# d is the input data, npls is the number of pls components, ntime is the
# number of time trends, beta0cov, and beta_x_cov are the covariance models
# for the intercept and the time trends, respectively and must be either
# "exp" or "iid", nu_re specifies whether to include a random effect in the nugget,
# range, sill, nugget, and random are vectors of initial values for each respective
# parameter (all vectors must be the same length), and grad_h means something
# but I'm not sure what.

model_fit_params <- function(d, x=NULL, npls=2, ntime=2, beta0cov="exp",
  beta_x_cov="exp", nu_re=TRUE, 
  cov_beta_nugget=TRUE, 
  cov_nu_nugget=TRUE, 
  range=c(3), sill=c(-5), 
  nugget=c(-5),
  random=c(-5), grad_h=0.00001, dftrend=8 * 19) {

  if (!beta0cov %in% c("exp", "iid")) {
      stop(paste("Incorrect specification of the covariance form for beta0. \n",
          " Should be either 'exp' or 'iid'\n"))
  }
  if (!beta_x_cov %in% c("exp", "iid")) {
      stop(paste("Incorrect specification of the covariance form for beta0.",
          " \n Should be either 'exp' or 'iid'\n"))
  }

  if (length(npls) == 1) {
      npls <- rep(npls, ntime + 1)
  }

  cov_beta <- list(covf = c(beta0cov, rep(beta_x_cov, ntime)),
     nugget = c(cov_beta_nugget, rep(cov_beta_nugget, ntime)))
  cov_nu <- list(covf = "exp", nugget = cov_nu_nugget, random.effect = nu_re)

  # Find the indices for the aqs and fixed monitors that correspond to
  # the column order of the st_matrix.
  aqs_and_fixed_id <- unique(d$native_id.monitors[d$aqs_and_fixed_indices])
  st_matrix_cols <- colnames(d$st_matrix)
  st_matrix_aqs_and_fixed_indices <- which(st_matrix_cols %in% aqs_and_fixed_id)

  # Get Smoothed Time Trend using only the fixed and aqs locations.
  stdata_trendonly <- createSTdata(obs = d$st_matrix[ ,st_matrix_aqs_and_fixed_indices],
      covars = d$lambert.monitors[d$aqs_and_fixed_indices, ],
      n.basis = ntime,
      df = dftrend, transform.obs=log)

  noFixedRows <- rowSums(is.na(d$st_matrix[ ,st_matrix_aqs_and_fixed_indices]))==
    dim(d$st_matrix[ ,st_matrix_aqs_and_fixed_indices])[2]
  if (sum(noFixedRows) > 0) {
    stop(paste("There are not any fixed or AQS measurements for these dates: ",
      paste(rownames(d$st_matrix)[which(noFixedRows)], sep=","),
      ", but there are other measurements for these dates. This is not allowed.", sep=""))
  }

  # Create PLS components
  pls_out <- create_pls_comps(d$st_matrix, st_matrix_aqs_and_fixed_indices,
    d$covars.monitors, stdata_trendonly$trend, npls = npls, ntime = ntime)

  st_covar_matrix <- as.data.frame(cbind(
    pls_out$pls.covars,
    d$lambert.monitors,
    d$lon_lat,
    as.data.frame(d$type.monitors)))
  rownames(st_covar_matrix) <- as.character(d$native_id.monitors)


  st_covar_matrix <- rename(st_covar_matrix, c("d$type.monitors"="type"))

  # Create primary STdata object with all of the observations. In the intro
  # this is called mesa.data.

  st_data <- createSTdata(obs = d$st_matrix, covars = st_covar_matrix, transform.obs=log)

  st_data$trend <- stdata_trendonly$trend
  st_data$trend.fnc <- stdata_trendonly$trend.fnc

  # Create list of names ID'ing the covars with each time trend
  lur <- list()
  for (b in 0:ntime) {
    lur[[b + 1]] <- paste("b", b, "comp", 1:npls[b + 1], sep = "")
  }

  # Create primary STmodel object. In the intro this is called mesa.model.
  locations <- list(coords = c("lambert_x", "lambert_y"),
                    long.lat = c("longitude", "latitude"),
                    others = "type")
  st_model <- createSTmodel(STdata = st_data, LUR = lur, cov.beta = cov_beta,
      cov.nu = cov_nu, scale = FALSE, locations = locations)
  st_model$locations$type <- d$type.monitors

  if( is.null(x) ) {
    x.init <- create_init(st_model, range=range, sill=sill,
      nugget=nugget, random=random)
  } else {
    x.init <- x
  }

  # estimate the parameters
  st_estimate <- estimate_edited.STmodel(object = st_model,
                          x = x.init,
                          x.fixed = NULL,
                          type="p",
                          hessian.all=FALSE,
                          method="L-BFGS-B",
                          h = grad_h,
                          diff.type=1,
                          lower=-15,
                          upper=15,
                          control=list(trace=3, factr=1e7, lmm=100))

  out <- list(st_data = st_data, st_model = st_model, pls_out = pls_out,
    st_estimate = st_estimate)
  return (out)
}

create_init <- function(st_model, range=c(0,-1), sill=c(0,-1), nugget=c(0,-1),
  random=c(0,-1)) {

  parms <- loglikeSTnames(st_model, all=FALSE)
  num.parms <- length(parms)
  num.runs <- length(range)
  x_init <- matrix(data = 0, nrow = num.parms, ncol = num.runs)

  for (i in 1:num.runs) {
    x_init[grep('range', parms), i] <- range[i]
    x_init[grep('sill', parms), i] <- sill[i]
    x_init[grep('nugget', parms), i] <- nugget[i]
    x_init[grep('random', parms), i] <- random[i]
  }
  print("Initial values:")
  print(x_init)

  return (x_init)
}


estimate_edited.STmodel <- function(object, x, x.fixed=NULL, type="p",
                             h=1e-3, diff.type=1, hessian.all=FALSE,
                             lower=-15, upper=15, method="L-BFGS-B",
                             control=list(trace=3, maxit=1000), ...){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")

  ##get size of the models
  dimensions <- loglikeSTdim(object)

  ##ensure lower case
  type <- tolower(type)
  ##first check that type is valid
  SpatioTemporal:::stCheckType(type)

  ##Second check the input x
  x <- as.matrix(x)
  ##if x has length one we need to construct a matrix of initial values
  if( length(x)==1 ){
    x <- matrix(seq(-5,5,length.out=x), dimensions$nparam.cov, x, byrow=TRUE)
  }
  ##check that starting point is valid
  tmp <- SpatioTemporal:::stCheckX(x, x.fixed, dimensions, type, object)
  x.all <- tmp$x.all
  x <- tmp$x
  x.fixed <- tmp$x.fixed

  ##Default values for control
  control <- defaultList(control, eval(formals(estimate_edited.STmodel)$control) )

  ##define local version of gradient function, fixing h and diff.type
  loglikeSTGrad.loc <- function(x, STmodel, type, x.fixed){
    loglikeSTGrad(x, STmodel, type, x.fixed, h=h, diff.type=diff.type)
  }##function loglikeSTGrad.loc

  ##attempt to fit the model for each of the provided starting values.
  res <- as.list( rep(NA, dim(x)[2]) )
  ##vector with convergence information and optimal values
  conv <- rep(FALSE, dim(x)[2])
  value <- rep(NA, dim(x)[2])

  ##make sure that likelihood evaluates
  err <- tryCatch( loglikeST(x[,1], STmodel=object, type=type,
                             x.fixed=x.fixed), silent=TRUE )
  if( inherits(err,"try-error") ){
    stop( paste("log-likelihood fails at first starting point with:\n",
                err[[1]]) )
  }

  #############
  # J. Keller edit; not in original package code
  #
  # Define local version of the likelihood function, to use in calculating
  # Hessians using numDeriv package
  #
  loglikeST.loc <- function(x) loglikeST(x, object)
  ############

  ##ensure that we are doing maximisation
  control$fnscale <- -1
  ##loop over starting values
  for(i in 1:dim(x)[2]){
    if( control$trace!=0 ){
      message( paste("Optimisation using starting value ",
                     i, "/", dim(x)[2], sep="") )
    }
    possibleError <- tryCatch(
      resTmp <- run_model(x[,i], x.all[,i], loglikeST, loglikeST.loc,
        loglikeSTGrad.loc, object, type, x.fixed,
        method, control, lower, upper),
      error=function(e) print(e)
    )
    if (inherits(possibleError, "error")) {
      if (exists("resTmp")) {
        if( control$trace!=0 ){
          message("Solver error. Starting again with final params + noise.")
        }
        resTmp <- run_model(jitter(resTmp$res$par), x.all[,i], loglikeST, loglikeST.loc,
          loglikeSTGrad.loc, object, type, x.fixed,
          method, control, lower, upper)
        } else {
          if( control$trace!=0 ){
            message("Solver error. Starting again with initial params + noise.")
          }
          resTmp <- run_model(jitter(x[,i]), x.all[,i], loglikeST, loglikeST.loc,
            loglikeSTGrad.loc, object, type, x.fixed,
            method, control, lower, upper)
        }
    }

    if (!exists("resTmp")) {
      if (resTmp$conv == FALSE) {
        if( control$trace!=0 ){
          message("Didn't converge. Starting again with final params + noise.")
        }
        resTmp <- run_model(jitter(resTmp$res$par), x.all[,i], loglikeST, loglikeST.loc,
          loglikeSTGrad.loc, object, type, x.fixed,
          method, control, lower, upper)
      }
  }

    res[[i]] <- resTmp$res
    conv[[i]] <- resTmp$conv
    value[[i]] <- resTmp$value

    #if (resTmp$conv == TRUE) {
    #  message("Converged. Skipping the rest of the starting points.")
    #  break
    #} else {
    #  message("Didn't converge again. Continuing to the next starting point.")
    #}
  }##for(i in 1:dim(x)[2])

  if( all(is.na(res)) ){
    stop("All optimisations failed, consider trying different starting values.")
  }

  ##extract summaries of the optimisations
  status <- data.frame(value=value, convergence=logical(length(res)),
                       conv=(conv==1))
  par.cov <- matrix(NA, dimensions$nparam.cov, length(res))
  par.all <- matrix(NA, dimensions$nparam, length(res))
  for(i in 1:length(res)){
    if( all(!is.na(res[[i]])) ){
      status$convergence[i] <- res[[i]]$convergence==0
      par.cov[,i] <- res[[i]]$par.cov$par
      par.all[,i] <- res[[i]]$par.all$par
    }
  }
  ##add names to the summaries
  rownames(par.all) <- loglikeSTnames(object, all=TRUE)
  rownames(par.cov) <- loglikeSTnames(object, all=FALSE)

  ##pick out the converged option with the best value
  Ind.overall <- which.max(value)
  if(any(conv==TRUE)){
    ##mark no-converged values as NA to avoid picking these
    value[!conv] <- NA
  }
  ##extract the best value
  Ind <- which.max(value)
  res.best <- res[[Ind]]

  ##collect status results
  summary <- list(status=status, par.all=par.all, par.cov=par.cov, x.fixed=x.fixed)

  if(hessian.all==TRUE){
    if(type!="f"){
      x.fixed <- res.best$par.all$fixed
      x <- res.best$par.all$par[ is.na(x.fixed) ]
      res.best$hessian.all <- loglikeSTHessian(x, object, type="f",
                                               x.fixed=x.fixed, h=h)
      ##standard error
      suppressWarnings( par.sd <- sqrt(-diag(solve(res.best$hessian.all))) )
      res.best$par.all$sd <- NA
      res.best$par.all$sd[ is.na(x.fixed) ] <- par.sd
      ##update t-statistic for the best result, all parameters
      res.best$par.all$tstat <- res.best$par.all$par/res.best$par.all$sd
    }else{
      ##replicate hessian for all parameters so output is consistent
      res.best$hessian.all <- res.best$hessian
    }
  }##if(hessian.all==TRUE)

  ##return result
  out <- list(res.best=res.best, res.all=res, summary=summary)
  class(out) <- "estimateSTmodel"

  return( out )
}##function estimate_edited.STmodel


run_model <- function(x, x.all, loglikeST, loglikeST.loc, loglikeSTGrad.loc,
  object, type, x.fixed,
  method, control, lower, upper) {

  conv <- FALSE
  value <- NA

  ##get size of the models
  dimensions <- loglikeSTdim(object)

  res <- optim(x, loglikeST, gr=loglikeSTGrad.loc,
            STmodel=object, type=type, x.fixed=x.fixed,
            method=method, control=control, hessian=TRUE,
            lower=lower, upper=upper)
  ##has optimisation converged?
  if( all( !is.na(res) ) ){
    ##then compute convergence criteria

    #############
  # J. Keller edit; not in original package code
  #
  # Calculate Hessian using numDeriv package
    res$optim_hessian <- res$hessian
    res$hessian <- hessian(loglikeST.loc, res$par) #, method.args=list(d=0.01, r=6))

    ###############

    conv <- (res$convergence==0 &&
                all(eigen(res$hessian)$value < -1e-10))
    ##extract ML-value
    value <- res$value

    ##add convergence and initial parameters
    res$conv <- conv
    res$par.cov <- data.frame(par=double(dimensions$nparam.cov), sd=NA,
                                   fixed=double(dimensions$nparam.cov),
                                   init=double(dimensions$nparam.cov),
                                   tstat=double(dimensions$nparam.cov))
    res$par.all <- data.frame(par=double(dimensions$nparam), sd=NA,
                                   fixed=double(dimensions$nparam),
                                   init=double(dimensions$nparam),
                                   tstat=double(dimensions$nparam))

    ##add standard deviations
    #
    # Edit by J. Kelller, January 2013
    # Original code had:
    #suppressWarnings( par.sd <- sqrt(-diag(solve(res$hessian))) )
    #
    # If optimization resulted in (computationally) singular hessian, then this
    # throws an error that is not picked up and execution halts.
    #Switched to this code:
    par.sd <- sqrt(-diag(tryCatch(solve(res$hessian), error=function(e) rep(-10000, length(res$par)))))


    ##initial value
    if( type!="f" ){
      par.type <- "par.cov"
    }else{
      par.type <- "par.all"
    }
    ##parameters
    res[[par.type]]$init <- x.all
    res[[par.type]]$par <- res[[par.type]]$fixed <- x.fixed
    res[[par.type]]$par[is.na(x.fixed)] <- res$par
    ##standard error
    res[[par.type]]$sd[is.na(x.fixed)] <- par.sd

    if( type!="f" ){
      ##compute regression parameters
      tmp <- predict(object, res$par.cov$par, only.pars=TRUE,
                     pred.var=FALSE, type=type)$pars
      res$par.all$par <- c(tmp$gamma.E, tmp$alpha.E, res$par.cov$par)
      N.reg <- length(tmp$gamma.E)+length(tmp$alpha.E)
      res$par.all$sd <- c(rep(NA,N.reg), res$par.cov$sd)
      res$par.all$init <- c(rep(NA,N.reg), x.all)
      res$par.all$fixed <- c(rep(NA,N.reg), x.fixed)
    }else{
      ##all the covariance parameters
      I <- (dimensions$nparam-dimensions$nparam.cov+1):dimensions$nparam
      res$par.cov <- res$par.all[I,,drop=FALSE]
    }
    ##t-statistic
    res$par.cov$tstat <- res$par.cov$par / res$par.cov$sd
    res$par.all$tstat <- res$par.all$par / res$par.all$sd

    ##add names to the variables.
    rownames(res$par.all) <- loglikeSTnames(object, all=TRUE)
    rownames(res$par.cov) <- loglikeSTnames(object, all=FALSE)

    if( type!="f" ){
      names(res$par) <- loglikeSTnames(object, all=FALSE)[is.na(x.fixed)]
    }else{
      names(res$par) <- loglikeSTnames(object, all=TRUE)[is.na(x.fixed)]
    }
  }##if(all(!is.na(res)))
  out <- list("res" = res, "conv" = conv, "value" = value)
  return(out)
}


###################################
## Functions for crossvalidation ##
###################################
##Functions in this file:
## estimateCV.STmodel          EX:ok
## estimateCV                  EX:S3 method
## print.estCVSTmodel          EX:ok
## summary.estCVSTmodel        EX:ok
## print.summary.estCVSTmodel  EX:with summary.estCVSTmodel
## coef.estCVSTmodel           EX:ok
## boxplot.estCVSTmodel        EX:ok

##' Functions that perform cross-validated parameter estimation and prediction
##' for the spatio-temporal model.
##'
##' For \code{predictCV.STmodel} the parameters used to compute predictions for the left
##' out observations can be either a single vector or a matrix.
##' For a single vector the same parameter values will be used for all
##' cross-validation predictions; for a matrix the parameters in \code{x[,i]}
##' will be used for the predictions of the i:th cross-validation set (i.e. for
##' \code{Ind.cv[,i]}). Suitable matrices are provided in the output from
##' \code{estimateCV.STmodel}.
##'
##' The cross-validation groups are given by \code{Ind.cv}. \code{Ind.cv} should
##' be either a (number of observations) - by - (groups) logical matrix or an
##' \emph{integer valued} vector with length equal to (number of observations).
##' If a matrix then each column defines a cross-validation set with the
##' \code{TRUE} values marking the observations to be left out. If a vector then
##' \code{1}:s denote observations to be dropped in the first cross-validation
##' set, \code{2}:s observations to be dropped in the second set, etc.
##' Observations marked by values \code{<=0} are never dropped. See
##' \code{\link{createCV}} for details.
##'
##' @title Cross-Validated Estimation and Prediction
##'
##' @param object \code{STmodel} object for which to perform cross-validation.
##' @param x Either a vector or matrix of starting point(s) for the optimisation,
##'   see \code{\link{estimate.STmodel}}; or a matrix with parameters, the i:th
##'   row being used for prediction of the i:th cross-validation set. For
##'   prediction either a \code{estCVSTmodel} or \code{estimateSTmodel} object,
##'   results from \code{\link{estimateCV.STmodel}} or
##'   \code{\link{estimate.STmodel}}, can be used.
##' @param Ind.cv \code{Ind.cv} defines the cross-validation scheme.  Either a
##'   (number or observations) - by - (groups) logical matrix or an \emph{integer
##'   valued} vector with length equal to (number or observations). For
##'   \code{predictCV.STmodel} \code{Ind.cv} can be infered from \code{x} if
##'   \code{x} is a \code{estCVSTmodel} object
##'   See further \code{\link{createCV}}.
##' @param control A list of control parameters for the optimisation.
##'   See \code{\link[stats:optim]{optim}} for details; setting \code{trace}=0
##'   eliminates all ouput.
##' @param verbose.res A \code{TRUE}/\code{FALSE} variable indicating if full
##'   results from \code{\link{estimate.STmodel}} for each CV group should be
##'   returned; defaults to \code{FALSE}
##' @param ... All additional parameters for \code{\link{estimate.STmodel}}
##'   or \code{\link{predict.STmodel}}.
##'   For \code{\link{predict.STmodel}} a number of parameters are set in
##'   \code{predictCV.STmodel} and can \strong{NOT} be overriden, these are
##'   \code{nugget.unobs}, \code{only.pars=FALSE}, and
##'   \code{combine.data=FALSE}.
##'
##' @return Either a \code{estCVSTmodel} object with elements:
##'   \item{status}{Data.frame with convergence information and best function
##'                 value for each cross-validation group.}
##'   \item{Ind.cv}{The cross-validation grouping.}
##'   \item{x.fixed}{Fixed parameters in the estimation, see
##'                  \code{\link{estimate.STmodel}}.}
##'   \item{x.init}{Matrix of inital values used, i.e. \code{x} from the input.}
##'   \item{par.all, par.cov}{Matrices with estimated parameters for each
##'                           cross-validation group.}
##'   \item{par.all.sd, par.cov.sd}{Standard deviations computed from the
##'     Hessian/information matrix for set of estimated parameters.}
##'   \item{res.all}{Estimation results for each cross-validation group,
##'                  contains the output from the \code{\link{estimate.STmodel}}
##'                  calls, only included if \code{verbose.res=TRUE}.}
##' Or a \code{predCVSTmodel} object with elements:
##'   \item{opts}{Copy of the \code{opts} field in the output from
##'               \code{\link{predict.STmodel}}.}
##'   \item{Ind.cv}{The cross-validation grouping.}
##'   \item{pred.obs}{A data.frame with a copy of observations from
##'                   \code{object$obs}, predictions (for different model
##'                   components), variances, and residuals. Variance field will
##'                   be missing if \code{pred.var=FALSE}.}
##'   \item{pred.all}{A list with time-by-location data.frames containing
##'                   predictions and variances for all space-time locations as
##'                   well as predictions and variances for the
##'                   beta-fields. Unobserved points are \code{NA} for the
##'                   option \code{only.obs=TRUE}.}
##'
##' @example Rd_examples/Ex_estimateCV_STmodel.R
##'
##' @author Johan Lindstr�m
##' @family STmodel methods
##' @family cross-validation functions
##' @family estCVSTmodel methods
##' @method estimateCV STmodel
##' @export
# If fullCV is TRUE, PLS and time trends will be re-fit for each member of the CV,
# using the specified hyperparameters. Otherwise, the hyperparameters will not be used.
estimateCV_edited.STmodel <- function(object, x, Ind.cv, control=list(trace=3),
                               verbose.res=FALSE,
                               fullCV=FALSE,
                               preprocessedData=NULL,
                               npls=2, 
                               ntime=2, 
                               beta0cov="exp",
                              beta_x_cov="exp",
                              nu_re=TRUE, 
                              cov_beta_nugget=TRUE, 
                              cov_nu_nugget=TRUE,
                              dftrend=8 * 19, ...){

  ##check class belonging
  stCheckClass(object, "STmodel", name="object")
  ##check cross-validation groups
  Ind.cv <- stCheckInternalCV(Ind.cv)

  ##Default values for control
  control <- defaultList(control, eval(formals(estimate_edited.STmodel)$control) )

  ##ensure that Ind.cv is a matrix
  Ind.cv <- as.matrix(Ind.cv)
  if( dim(Ind.cv)[2]==1 ){
    N.CV.sets <- max(Ind.cv,na.rm=TRUE)
  }else{
    N.CV.sets <- dim(Ind.cv)[2]
  }

  ##get size of the models
  dimensions <- loglikeSTdim(object)

  res <- list()
  for(i in 1:N.CV.sets){
    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==i
    }else{
      Ind.current <- as.logical( Ind.cv[,i] )
    }

#--> time trends for object.aux
#--> PLS for object.aux

    ##lets estimate parameters for this set
    if( control$trace!=0 ){
      message( "\n***************************")
      message( paste("Estimation of CV-set ", i, "/", N.CV.sets, sep="") )
    }


    ############ C Tessum edit
    if( fullCV==TRUE ) {
      pd_current <- preprocessedData
      # Drop locations
      drop_IDs <- levels(as.factor(object$obs$ID[which(Ind.current)]))
      keep_idx <- which(!colnames(pd_current$st_matrix) %in% drop_IDs)
      pd_current$st_matrix <- pd_current$st_matrix[ , keep_idx]
      pd_current$covars.monitors <- pd_current$covars.monitors[keep_idx, ]
      pd_current$lon_lat <- pd_current$lon_lat[keep_idx, ]
      pd_current$lambert.monitors <- pd_current$lambert.monitors[keep_idx, ]
      pd_current$type.monitors <- pd_current$type.monitors[keep_idx]
      pd_current$native_id.monitors <- pd_current$native_id.monitors[keep_idx]
      pd_current$aqs_and_fixed_indices <- which(pd_current$type.monitors == "AQS" | pd_current$type.monitors == "FIXED")

      # Run the model.
       res[[i]] <- model_fit_params(pd_current, x=x, npls=npls, ntime=ntime, beta0cov=beta0cov,
          beta_x_cov=beta_x_cov, nu_re=nu_re, cov_beta_nugget=cov_beta_nugget,
          cov_nu_nugget=cov_nu_nugget, dftrend=dftrend, ...)$st_estimate
    } else {
      ##create data matrices that contains observations
      object.aux <- dropObservations(object, Ind.current)
      res[[i]] <- estimate_edited.STmodel(object=object.aux, x=x, control=control, ...)
    }
    ############ End C Tessum edit
  }##for(i in 1:dim(Ind.cv)[2])

  ##status of optimisations
  status <- data.frame(value=sapply(res, function(x){x$res.best$value}),
                       convergence=sapply(res, function(x){x$res.best$convergence==0}),
                       conv=sapply(res, function(x){x$res.best$conv}))
  ##add hessian eigen-values to status
  tmp <- sapply(res, function(x){ range( -eigen(x$res.best$hessian)$value) })
  status$eigen.min <- tmp[1,]
  if( is.null( res[[1]]$res.best$hessian.all) ){
    status$eigen.all.min <- NA
  }else{
    tmp <- sapply(res,
                  function(x){ range( -eigen(x$res.best$hessian.all)$value) })
    status$eigen.all.min <- tmp[1,]
  }

  ##matrices with the estimates parameters (accounting for x.fixed)
  par.cov <- matrix(NA, dim(res[[1]]$res.best$par.cov)[1], N.CV.sets)
  par.cov.sd <- matrix(NA, dim(res[[1]]$res.best$par.cov)[1], N.CV.sets)
  par.all <- matrix(NA, dim(res[[1]]$res.best$par.all)[1], N.CV.sets)
  par.all.sd <- matrix(NA, dim(res[[1]]$res.best$par.all)[1], N.CV.sets)

  ##best estimates for each CV-group
  for(i in 1:N.CV.sets){
    par.cov[,i] <- res[[i]]$res.best$par.cov$par
    par.cov.sd[,i] <- res[[i]]$res.best$par.cov$sd
    par.all[,i] <- res[[i]]$res.best$par.all$par
    par.all.sd[,i] <- res[[i]]$res.best$par.all$sd
  }
  ##all initial values
  x.init <- sapply(res[[1]]$res.all, function(x){ x$par.all$init})
  ##also return a list with all optimisation results?
  if( verbose.res ){
    res.all <- res
  }else{
    res.all <- NULL
  }

  rownames(par.cov) <- rownames(res[[1]]$res.best$par.cov)
  rownames(par.cov.sd) <- rownames(res[[1]]$res.best$par.cov)
  rownames(par.all) <- rownames(res[[1]]$res.best$par.all)
  rownames(par.all.sd) <- rownames(res[[1]]$res.best$par.all)
  rownames(x.init) <- rownames(res[[1]]$res.best$par.all)
  ##drop NA's from x.init (occurs if only covariance parameters where given)
  x.init <- x.init[!apply(is.na(x.init),1,all),,drop=FALSE]
  ##Return the estimated parameters from the cross-validation
  out <- list(par.cov=par.cov, par.cov.sd=par.cov.sd,
              par.all=par.all, par.all.sd=par.all.sd,
              res.all=res.all, status=status, Ind.cv=Ind.cv,
              x.fixed=res[[1]]$summary$x.fixed,
              x.init=x.init)
  class(out) <- "estCVSTmodel"

  return( out )
}##function estimateCV_edited.STmodel

stCheckInternalCV <- function(Ind.cv, Icv.vector=TRUE){
  if( is.null(Ind.cv) ){
    stop("Ind.cv must be defined (is NULL)")
  }else if( is.vector(Ind.cv) ){
    ##ensure integers
    Ind.cv <- as.integer(Ind.cv)
    uInd.cv <- sort(unique(Ind.cv))
    ##find missing/empty groups
    I.rep <- match(1:max(uInd.cv),uInd.cv)
    I.drop <- which(is.na(I.rep))
    if( length(I.drop!=0) ){
      ##renumber CV groups
      I.rep[ !is.na(I.rep) ] <- 1:sum(!is.na(I.rep))
      I.rep <- c(I.rep,0)
      Ind.cv[Ind.cv==0] <- length(I.rep)
      Ind.cv <- I.rep[Ind.cv]
    }
  }else{
    I.drop <- which( colSums(Ind.cv)==0 )
    ##renumber CV groups
    if( length(I.drop!=0) ){
      Ind.cv <- Ind.cv[, -I.drop, drop=FALSE]
    }
    if( Icv.vector && max(rowSums(Ind.cv))==1 ){
      Ind.cv <- apply(Ind.cv, 1, function(x){ x=which(x); if(length(x)!=1){x=0};
                                              return(x)})
    }
  }##if( is.vector(Ind.cv) ){...}else{...}
  if( length(I.drop!=0) ){
    warning( paste("Empty cross validation group(s) dropped:",
                   paste(I.drop, collapse=", ")) )
  }
  return(Ind.cv)
}##function stCheckInternalCV
