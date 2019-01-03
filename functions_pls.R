# libraries
library(pls)

pls_regression <- function(covar.matrix, ys, npls) {
  pls.data <- c()
  mask <- !is.na(covar.matrix)
  covar.matrix <- ifelse(mask, covar.matrix, 0)
  pls.data$X <- covar.matrix
  pls.data$y <- ys
  pls.obj <- plsr(y ~ X, ncomp = npls, validation = "none", data = pls.data, na.action = na.fail)
  return(pls.obj)
}

pls_predict <- function(pls.obj, covar.matrix) {
  targets.data <- c()
  mask <- !is.na(covar.matrix)
  covar.matrix <- ifelse(mask, covar.matrix, 0)
  targets.data$X <- covar.matrix
  targets.scores <- predict(pls.obj, newdata = targets.data, type = "scores")
  return(targets.scores)
}



# Creates the PLS components
create_pls_comps <- function(st_matrix, aqs_and_fixed_indices, covars.monitors, trend,
  npls, ntime) {
  # Number of time trends
  if (dim(st_matrix)[1] != dim(trend)[1]) {
    stop(paste("create_pls_comps: st_matrix has a different number of rows than trend: ",
      dim(st_matrix)[1], "!=", dim(trend)[1], sep = ""))
  }
  ntrend <- ncol(trend) - 1
  all.fixed <- colnames(st_matrix)[aqs_and_fixed_indices]
  # Fixed Sites Regress fixed sites' time series on the trend functions To generate
  # 'empirical betas'
  emp.beta.matrix <- matrix(NA, length(all.fixed), ntrend + 1)
  rownames(emp.beta.matrix) <- all.fixed
  for (j in all.fixed) {
    jmiss <- is.na(st_matrix[, j])
    if (length(st_matrix[!jmiss, j]) <= max(npls)) {
      stop(paste("in create_pls_comps: number of measurements <= npls for monitor ID ",
        colnames(st_matrix)[j], ", can't make a linear regression", sep = ""))
    }
    temp <- lm(st_matrix[!jmiss, j] ~ as.matrix(trend[!jmiss, 1:ntime]), na.action = na.omit)
    emp.beta.matrix[j, ] <- coef(temp)
  }
  if (anyNA(emp.beta.matrix)) {
    stop("in pls_regression: emp.beta.matrix contains an NA value")
  }
  # Do PLS regression of empirical betas on the GIS covariates
  covar.matrix <- scale(as.matrix(covars.monitors[all.fixed, ]))

  #
  pls.out <- list()
  for (j in 1:(ntime + 1)) {
    pls.out[[j]] <- pls_regression(covar.matrix, emp.beta.matrix[, j], npls[j])
  }

  all.nonfixed <- setdiff(colnames(st_matrix), all.fixed)
  nonfixed.pls.out <- list()
  if (length(all.nonfixed) >= 1) {
    nonfixed.covars <- as.matrix(covars.monitors[all.nonfixed, ])
    nonfixed.covars <- scale(nonfixed.covars, scale = attr(covar.matrix, "scaled:scale"),
      center = attr(covar.matrix, "scaled:center"))

    for (j in 1:(ntime + 1)) {
      nonfixed.pls.out[[j]] <- pls_predict(pls.out[[j]], nonfixed.covars)
    }
    # construct the PLS covariate matrix
    pls.covars <- rbind(pls.out[[1]]$scores, nonfixed.pls.out[[1]])
  } else {
    pls.covars <- pls.out[[1]]$scores
  }

  rownames(pls.covars) <- c(all.fixed, all.nonfixed)
  colnames(pls.covars)[1:npls[1]] <- paste("b0comp", 1:npls[1], sep = "")
  for (j in 2:(ntime + 1)) {
    if (length(all.nonfixed) != 0) {
      pls.covars <- cbind(pls.covars, rbind(pls.out[[j]]$scores, nonfixed.pls.out[[j]]))
    } else {
      pls.covars <- cbind(pls.covars, pls.out[[j]]$scores)
    }
    tot.cols <- ncol(pls.covars)
    colnames(pls.covars)[(tot.cols - npls[j] + 1):tot.cols] <- paste("b", j -
      1, "comp", 1:npls[j], sep = "")
  }

  if (anyNA(pls.covars)) {
    stop("create_pls_comps: pls.covars contains an NA value")
  }

  # Reorder covariates so IDs match st_matrix
  pls.covars <- pls.covars[match(colnames(st_matrix),rownames(pls.covars)), ]

  return(list(pls.out = pls.out, pls.covars = pls.covars, covar.scale = attr(covar.matrix,
    "scaled:scale"), covar.center = attr(covar.matrix, "scaled:center")))
}
