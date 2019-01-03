
runCV <- function(st_model, st_estimate, include_comco=TRUE, include_home=TRUE) {

  out <- list()

  x.init <- coef(st_estimate, pars="cov")[,c("par","init")]

  Ind.cv.fixed <- createCV(st_model, groups=10, min.dist=.1, option="fixed",Icv.vector = F)
  st_cv_est_fixed <- estimateCV_edited.STmodel(st_model, x.init, Ind.cv.fixed)
  out$st_cv_est_fixed <- st_cv_est_fixed
  out$st_cv_predict_fixed <- predictCV(st_model, st_cv_est_fixed, LTA=TRUE, silent=FALSE)
  out$st_cv_predict_fixed_log <- predictCV(st_model, st_cv_est_fixed, LTA=TRUE, silent=FALSE,transform="unbiased")

  if (include_comco==TRUE) {
    comco_subset <- st_model$locations$ID[which(st_model$locations$type %in% c("COMCO"))]
    Ind.cv.comco <- createCV(st_model, groups=10, min.dist=.1, option="all", subset=comco_subset, Icv.vector = F)
    st_cv_est_comco <- estimateCV_edited.STmodel(st_model, x.init, Ind.cv.comco)
    out$st_cv_est_comco <- st_cv_est_comco
    out$st_cv_predict_comco <- predictCV(st_model, st_cv_est_comco, LTA=TRUE, silent=FALSE)
    out$st_cv_predict_comco_log <- predictCV(st_model, st_cv_est_comco, LTA=TRUE, silent=FALSE,transform="unbiased")
  }

  if (include_comco==TRUE) {
    Ind.cv.all <- createCV(st_model, groups=10, min.dist=.1, option=c("all"), Icv.vector = F)
    st_cv_est_all <- estimateCV_edited.STmodel(st_model, x.init, Ind.cv.all)
    out$st_cv_est_all <- st_cv_est_all
    out$st_cv_predict_all <- predictCV(st_model, st_cv_est_all, LTA=TRUE, silent=FALSE)
    out$st_cv_predict_all_log <- predictCV(st_model, st_cv_est_all, LTA=TRUE, silent=FALSE,transform="unbiased")
  }

  return(out)
  }
