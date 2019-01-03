# st_plots is a function that makes plots of a SpatioTemporal model and
# results. fname is the name of the saved data object to load,
# ids is a list of 2 fixed or AQS monitor locations to plot,
# homeID is the ID of one COMCO monitor to plot,
# covar is the name of a covariate to plot, and
# includeComco indicates whether to include COMCO plots.
# if includeComco is false, comcoID is ignored.
# includeHome indicates whether to include HOME plots.
# if includeHome is false, homeID is ignored.
st_plots <- function(fname, ids, comcoID, homeID, covar, includeComco, includeHome) {

load(fname)

print("fig2")
# https://cran.r-project.org/web/packages/SpatioTemporal/vignettes/ST_intro.pdf
# Figure 2: Counterclockwise from the top: Space-time locations of our observations divided
# into AQS (black) and MESA fixed (red) locations, qqnorm-plot for all observations, and
# dependence between observations and distance to coast.
png(filename=paste(fname,".fig2.png", sep=""))
layout(matrix(c(1,2,1,3), 2, 2))
par(mar=c(2.3,3.3,2,1), mgp=c(2,1,0))
plot(modelResults$st_data, "loc", main="Occurrence of Observations", xlab="", ylab="Location", col=c("black", "red"), legend.loc=NULL)
par(mar=c(3.3,3.3,2,1))
qqnorm(modelResults$st_data, line=1)
scatterplot(modelResults$st_data, covar=covar, xlab=covar, ylab="NOx (log ppb)", pch=19, cex=.25, smooth.args=list(span=4/5,degree=2))
dev.off()

print("fig5")
# https://cran.r-project.org/web/packages/SpatioTemporal/vignettes/ST_intro.pdf
# Figure 5: The smooth temporal trends fitted to data at one site. Fitted trends, residuals,
# auto-correlation, and partial auto-correlation functions are shown.
png(filename=paste(fname,".fig5.", ids[1], ".png", sep=""))
par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow=TRUE))
plot(modelResults$st_data, "obs", ID=ids[1], xlab="", ylab="NOx (log ppb)", main=paste("Temporal trend ",ids[1], sep=""))
plot(modelResults$st_data, "res", ID=ids[1], xlab="", ylab="NOx (log ppb)")
plot(modelResults$st_data, "acf", ID=ids[1])
# plot(modelResults$st_data, "pacf", ID=ids[1]) # This one doesn't work
dev.off()

print("fig6")
# https://cran.r-project.org/web/packages/SpatioTemporal/vignettes/ST_intro.pdf
# Figure 6: Comparing the fit of 2 smooth (black) and 3 deterministic (red) basis functions at
# two locations.
png(filename=paste(fname,".fig6.", ids, ".png", sep=""))
mesa.data.fnc <- updateTrend(modelResults$st_data, fnc=function(x){
  x = 2*pi*as.numeric(x)/365;
  return( cbind(x, sin(x), cos(x)) )
})
par(mfrow=c(2,1), mar=c(2.3,3.3,1.5,1), mgp=c(2,1,0))
for(i in ids) {
  plot(modelResults$st_data, ID=i, pch=c(19,NA), cex=.25, xlab="", ylab="NOx (log ppb)", main=paste("AQS site",i))
  plot(mesa.data.fnc, ID=i, add=TRUE, col=2, pch=NA)
}
dev.off()


print("fig11 fixed ")
# Figure 11: The two top panes shows predictions for a left-out location, in log and original
# scale (cf. Figure 1). In the left pane of the third row predictions have been plotted against
# observations; the points are coloured by location and grouping of data from single locations can
# be seen. On the right predictions and observations of temporal averages (35) are compared.
# In the last row, the left pane shows a QQ-plot for normalised residuals for the Gaussian
# model, coloured by season; the solid line gives the theoretical behaviour of N (0, 1) residuals.
# On the right, residuals (coloured by season) are plotted as a function of distance to A1-
# roads; smooths, for each season and all data, have been added to help identify any remaining
# structure.
png(filename=paste(fname,".fig11.", ids[1], ".fixed.png", sep=""))
par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
layout(matrix(c(1,1,2,2,3,4,5,6), 4, 2, byrow=TRUE))
plot(cvResults$st_cv_predict_fixed, ID=ids[1], xlab="", ylab="NOx (log ppb)", main=paste("Predictions for ", ids[1], sep=""), lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
plot(cvResults$st_cv_predict_fixed, ID=ids[1], pred.type="EX.mu", lty=4, lwd=2, col="blue", add=TRUE)
plot(cvResults$st_cv_predict_fixed, ID=ids[1], pred.type="EX.mu.beta", lty=2, lwd=2, col="green", add=TRUE)
plot(cvResults$st_cv_predict_fixed_log, ID=ids[1], xlab="", ylab="NOx (ppb)", main=paste("Predictions for ", ids[1], sep=""), pred.type="EX.pred", lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
plot(cvResults$st_cv_predict_fixed_log, ID=ids[1], pred.type="EX.mu", lty=4, lwd=2, col="blue", add=TRUE)
plot(cvResults$st_cv_predict_fixed_log, ID=ids[1], pred.type="EX.mu.beta", lty=2, lwd=2, col="green", add=TRUE)
legend("topright", c("Observations", "Predictions", "Contribution from beta", "Contribution from mean", "95% CI"), bty="n",
  lty=c(NA,1,2,4,NA), lwd=c(NA,2,2,2,NA), pch=c(19,NA,NA,NA,15), pt.cex=c(.75,NA,NA,NA,2.5),
  col=c("red", "black", "green", "blue", "grey"))
plot(cvResults$st_cv_predict_fixed, "obs", ID="all", pch=c(19,NA), cex=.25, lty=c(NA,2),
  col=c("ID", "black", "grey"), xlab="Observations",
  ylab="Predictions", main="Cross-validation NOx (log ppb)")
with(cvResults$st_cv_predict_fixed_log$pred.LTA, plotCI(obs, EX.pred, uiw=1.96*sqrt(VX.pred), xlab="Observations", ylab="Predictions", main="Temporal average NOx (ppb)"))
abline(0, 1, col="grey")
I.season <- as.factor(as.POSIXlt(cvResults$st_cv_predict_fixed$pred.obs$date)$mon+1)
levels(I.season) <- c(rep("Winter",2), rep("Spring",3), rep("Summer",3), rep("Fall",3), "Winter")
qqnorm(cvResults$st_cv_predict_fixed, norm=TRUE, main="Normalised residuals", col=I.season)
legend("bottomright", legend=as.character(levels(I.season)), pch=1, col=1:nlevels(I.season), bty="n")
scatterPlot(cvResults$st_cv_predict_fixed, STdata=modelResults$st_model, covar=covar, group=I.season, col=c(2:5,1), type="res",
  xlab=covar, ylab="Residuals", main="Residuals (log ppb)")
dev.off()

if (includeComco==TRUE) {
  print("fig11 comco ")
  png(filename=paste(fname,".fig11.comco.png", sep=""))
  par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
  layout(matrix(c(1,1,2,2,3,4,5,6), 4, 2, byrow=TRUE))
  plot(cvResults$st_cv_predict_comco, ID=comcoID, xlab="", ylab="NOx (log ppb)", main=paste("Predictions for ", comcoID, sep=""), lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
  plot(cvResults$st_cv_predict_comco, ID=comcoID, pred.type="EX.mu", lty=4, lwd=2, col="blue", add=TRUE)
  plot(cvResults$st_cv_predict_comco, ID=comcoID, pred.type="EX.mu.beta", lty=2, lwd=2, col="green", add=TRUE)
  plot(cvResults$st_cv_predict_comco_log, ID=comcoID, xlab="", ylab="NOx (ppb)", main=paste("Predictions for ", comcoID, sep=""), pred.type="EX.pred", lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
  plot(cvResults$st_cv_predict_comco_log, ID=comcoID, pred.type="EX.mu", lty=4, lwd=2, col="blue", add=TRUE)
  plot(cvResults$st_cv_predict_comco_log, ID=comcoID, pred.type="EX.mu.beta", lty=2, lwd=2, col="green", add=TRUE)
  legend("topright", c("Observations", "Predictions", "Contribution from beta", "Contribution from mean", "95% CI"), bty="n",
    lty=c(NA,1,2,4,NA), lwd=c(NA,2,2,2,NA), pch=c(19,NA,NA,NA,15), pt.cex=c(.75,NA,NA,NA,2.5),
    col=c("red", "black", "green", "blue", "grey"))
  plot(cvResults$st_cv_predict_comco, "obs", ID="all", pch=c(19,NA), cex=.25, lty=c(NA,2),
    col=c("ID", "black", "grey"), xlab="Observations",
    ylab="Predictions", main="Cross-validation NOx (log ppb)")
  with(cvResults$st_cv_predict_comco_log$pred.LTA, plotCI(obs, EX.pred, uiw=1.96*sqrt(VX.pred), xlab="Observations", ylab="Predictions", main="Temporal average NOx (ppb)"))
  abline(0, 1, col="grey")
  I.season <- as.factor(as.POSIXlt(cvResults$st_cv_predict_comco$pred.obs$date)$mon+1)
  levels(I.season) <- c(rep("Winter",2), rep("Spring",3), rep("Summer",3), rep("Fall",3), "Winter")
  qqnorm(cvResults$st_cv_predict_comco, norm=TRUE, main="Normalised residuals", col=I.season)
  legend("bottomright", legend=as.character(levels(I.season)), pch=1, col=1:nlevels(I.season), bty="n")
  scatterPlot(cvResults$st_cv_predict_comco, STdata=modelResults$st_model, covar=covar, group=I.season, col=c(2:5,1), type="res",
    xlab=covar, ylab="Residuals", main="Residuals (log ppb)")
  dev.off()
}

if (includeHome==TRUE) {
  print("fig11 home ")
  png(filename=paste(fname,".fig11.home.png", sep=""))
  par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
  layout(matrix(c(1,1,2,2,3,4,5,6), 4, 2, byrow=TRUE))
  plot(cvResults$st_cv_predict_home, ID=homeID, xlab="", ylab="NOx (log ppb)", main=paste("Predictions for ", homeID, sep=""), lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
  plot(cvResults$st_cv_predict_home, ID=homeID, pred.type="EX.mu", lty=4, lwd=2, col="blue", add=TRUE)
  plot(cvResults$st_cv_predict_home, ID=homeID, pred.type="EX.mu.beta", lty=2, lwd=2, col="green", add=TRUE)
  plot(cvResults$st_cv_predict_home_log, ID=homeID, xlab="", ylab="NOx (ppb)", main=paste("Predictions for ", homeID, sep=""), pred.type="EX.pred", lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
  plot(cvResults$st_cv_predict_home_log, ID=homeID, pred.type="EX.mu", lty=4, lwd=2, col="blue", add=TRUE)
  plot(cvResults$st_cv_predict_home_log, ID=homeID, pred.type="EX.mu.beta", lty=2, lwd=2, col="green", add=TRUE)
  legend("topright", c("Observations", "Predictions", "Contribution from beta", "Contribution from mean", "95% CI"), bty="n",
    lty=c(NA,1,2,4,NA), lwd=c(NA,2,2,2,NA), pch=c(19,NA,NA,NA,15), pt.cex=c(.75,NA,NA,NA,2.5),
    col=c("red", "black", "green", "blue", "grey"))
  plot(cvResults$st_cv_predict_home, "obs", ID="all", pch=c(19,NA), cex=.25, lty=c(NA,2),
    col=c("ID", "black", "grey"), xlab="Observations",
    ylab="Predictions", main="Cross-validation NOx (log ppb)")
  with(cvResults$st_cv_predict_home_log$pred.LTA, plotCI(obs, EX.pred, uiw=1.96*sqrt(VX.pred), xlab="Observations", ylab="Predictions", main="Temporal average NOx (ppb)"))
  abline(0, 1, col="grey")
  I.season <- as.factor(as.POSIXlt(cvResults$st_cv_predict_home$pred.obs$date)$mon+1)
  levels(I.season) <- c(rep("Winter",2), rep("Spring",3), rep("Summer",3), rep("Fall",3), "Winter")
  qqnorm(cvResults$st_cv_predict_home, norm=TRUE, main="Normalised residuals", col=I.season)
  legend("bottomright", legend=as.character(levels(I.season)), pch=1, col=1:nlevels(I.season), bty="n")
  scatterPlot(cvResults$st_cv_predict_home, STdata=modelResults$st_model, covar=covar, group=I.season, col=c(2:5,1), type="res",
    xlab=covar, ylab="Residuals", main="Residuals (log ppb)")
  dev.off()
}

}
