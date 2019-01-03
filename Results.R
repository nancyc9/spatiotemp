source("C:\\Users\\Mei\\Dropbox\\wangmei\\MobileModel\\function_R2cv.R")

###Report AQS+fixed sites
# CV results, prediction log scale, pred. LTA for long term average, pred. obs for 2 week predictions 
# calculate R2 based on the observed and prediction 


obs.LTA.AF <- cvResults$st_cv_predict_fixed_log$pred.LTA$obs
pred.LTA.AF <- cvResults$st_cv_predict_fixed_log$pred.LTA$EX
pred.AF <- cvResults$st_cv_predict_fixed_log$pred.obs$EX
AF <- cvResults$st_cv_predict_fixed_log$pred.obs[complete.cases(pred.AF),]
obs.AF <- AF$obs
pred.AF <- AF$EX
R2cv_AF <- R2cv(obs.AF, pred.AF, log = F)
R2cv.LTA_AF <- R2cv(obs.LTA.AF, pred.LTA.AF, log = F)
R2reg_AF <- summary(lm(obs.AF~pred.AF))$r.squared
R2reg.LTA_AF <- summary(lm(obs.LTA.AF~pred.LTA.AF))$r.squared
# 
R2cv_AF
# 
R2cv.LTA_AF
# 
R2reg_AF
# 
R2reg.LTA_AF
# make a table for results of the Agency predictions 2 week Observations 
Agencypredictions <- matrix(c(2.35, 0.67, 0.77, 1.76, 0.71, 0.79 ),ncol=3,byrow=TRUE)
colnames(Agencypredictions) <- c("RMSE","R2 CV","R2 CV reg")
rownames(Agencypredictions) <- c("Agency","Agency + Shinyei")
Agencypredictions <- as.table(Agencypredictions) 
Agencypredictions

# make a table for results of the Agency predictions LTA 
Agencypredictions2 <- matrix(c(1.42, 0.27, 0.63, 1.35, 0.32, 0.48),ncol=3,byrow=TRUE)
colnames(Agencypredictions2) <- c("RMSE","R2 CV","R2 CV reg")
rownames(Agencypredictions2) <- c("Agency","Agency w/ Shinyei")
Agencypredictions2 <- as.table(Agencypredictions2) 
Agencypredictions2


### 

obs.LTA.comco <- cvResults$st_cv_predict_comco_log$pred.LTA$obs
pred.LTA.comco <- cvResults$st_cv_predict_comco_log$pred.LTA$EX
pred.comco <- cvResults$st_cv_predict_comco_log$pred.obs$EX
comco <- cvResults$st_cv_predict_comco_log$pred.obs[complete.cases(pred.comco),]
obs.comco <- comco$obs
pred.comco <- comco$EX
R2cv_comco <- R2cv(obs.comco, pred.comco, log = F)
R2cv.LTA_comco <- R2cv(obs.LTA.comco, pred.LTA.comco, log = F)
R2reg_comco <- summary(lm(obs.comco~pred.comco))$r.squared
R2reg.LTA_comco <- summary(lm(obs.LTA.comco~pred.LTA.comco))$r.squared
# 
R2cv_comco
# 
R2cv.LTA_comco
# 
R2reg_comco
# 
R2reg.LTA_comco


### 

### 

obs.LTA.All <- cvResults$st_cv_predict_All_log$pred.LTA$obs
pred.LTA.All <- cvResults$st_cv_predict_All_log$pred.LTA$EX
pred.All <- cvResults$st_cv_predict_All_log$pred.obs$EX
All <- cvResults$st_cv_predict_All_log$pred.obs[complete.cases(pred.All),]
obs.All <- All$obs
pred.All <- All$EX

R2cv_comco <- R2cv(obs.comco, pred.comco, log = F)


R2cv_All <- R2cv(obs.All, pred.All, log = F)
R2cv.LTA_All <- R2cv(obs.LTA.All, pred.LTA.All, log = F)
R2reg_all <- summary(lm(obs.All~pred.All))$r.squared
R2reg.LTA_All <- summary(lm(obs.LTA.All~pred.LTA.All))$r.squared
# 
R2cv_All
# 
R2cv.LTA_All
# 
R2reg_All
# 
R2reg.LTA_All





### 

require(ggplot2)

setwd("C:/Users/Mei/Desktop")
svg("Shinyei_AQS_2week_All.svg")
ggplot(AF,aes(obs,EX))+ geom_point(size=2)+ylim(0,40)+xlim(0,40)+geom_abline(slope=1,intercept=0)+coord_fixed()+xlab("Observation") + ylab("Prediction")+ggtitle("PM2.5 2-week Observations at Agency site: 10-fold CV on Agency+Shinyei")+  theme_classic()
dev.off()

svg("Shinyei_AQS_LTA_All.svg")
ggplot(cvResults$st_cv_predict_fixed_log$pred.LTA,aes(obs,EX))+ geom_point(size=3)+ylim(0,15)+xlim(0,15)+geom_abline(slope=1,intercept=0)+coord_fixed()+xlab("Observation") + ylab("Prediction")+ggtitle("PM2.5 Long-term Averages at Agency site: Agency+Shinyei")+  theme_classic()
dev.off()