#setwd("C:\\Users\\Mei\\Dropbox\\wangmei\\low_cost")
setwd("~/Desktop/MeiCode") ### wd for laptop 

load("dr0281new.RData")
### Agency cov
Agency<-read.csv("AQS_2wkavg.csv")
Agency$native_id <- as.factor(Agency$native_id)

##
### where is the agency_from_gis data? 
## 


Agency_gis <- agency_frm_gis


##remove portland site
Agency_gis<-Agency_gis[which(Agency_gis$latitude>46),]
a<-Agency$native_id
AQS_cov <- Agency_gis[Agency_gis$native_id %in% a, ]
dim(AQS_cov)
#18 831

### Fixed cov
Fixed_gis<-pscaa_gis
Fixed_gis<-Fixed_gis[which(Fixed_gis$latitude>46),]
Fixed<-read.csv("Fixed_2wkavg.csv")
Fixed$native_id <-as.factor(Fixed$native_id)
b<-Fixed$native_id
Fixed_cov <- Fixed_gis[Fixed_gis$native_id %in%  b, ]
dim(Fixed_cov)
#13 867

### Low_cost cov
Low_cost<-read.csv("low_cost_shinyei_2wkavg.csv")
Low_cost$native_id <-as.factor(Low_cost$native_id)
Low_cost_gis<-low_cost_gis
c<-Low_cost$native_id
Low_cost_cov <- Low_cost_gis[Low_cost_gis$native_id_new %in%  c, ]
Low_cost_cov$native_id <- Low_cost_cov$native_id_new
dim(Low_cost_cov)
#51 868

###Uniform the GIS cov
Low_cost_cov <- Low_cost_cov[,colnames(Low_cost_cov) %in%  colnames(AQS_cov) ]
Fixed_cov <- Fixed_cov[,colnames(Fixed_cov) %in%  colnames(AQS_cov) ]
AQS_cov<-AQS_cov[,colnames(AQS_cov) %in%  colnames(Low_cost_cov) ]

# load packages
#install.packages("SpatioTemporal")
#install.packages("numDeriv")
#install.packages("plyr")
#install.packages("plotrix")
#devtools::load_all("SpatioTemporal")
library("SpatioTemporal")
library("numDeriv")
library("plyr")
library("plotrix")

#install.packages("pls")
library(pls)

# Load functions
source("Data Preprocessing/functions_timeseries.R")
source("Data Preprocessing/func_covariate_preprocess.R")
source("Data Preprocessing/functions_misc_covar.R")
source("Data Preprocessing/func_covariate.R") 
source("Data Preprocessing/functions_model.R")
source("Data Preprocessing/functions_pls.R")
source("Data Preprocessing/function_cv_bak.R") 
source("Data Preprocessing/function_plots.R")


## parameters 
pollutant <- "PM25"
region <- "Seattle"
npls<-3 #3
ntime<-1
beta0cov<-"exp"
beta_x_cov<-"exp"
nu_re<-TRUE  #nugget random effect
dftrend<-4*19
cov_beta_nugget<-TRUE # can be false
cov_nu_nugget<-TRUE # can be true, false, or "type"
outputDir <- "output_Shinyei/"
include_comco <- TRUE

## Reading in the arguments
args <- commandArgs(TRUE)

print(args)
## Parse the arguments
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for (i in 1:length(args)){
     eval(parse(text=args[i]))
  }
}


if (!exists("pollutant")){
	stop(" Pollutant not specified. \n You must specify region using the arugment 'pollutant= '. \n options: NOx, NO2, PM25.")
}
if (!exists("ntime")){
	cat(" Number of time trends not specified. \n Use argument 'ntime= ' to specify the number of time trends. \n Defaulting to 1 time trend.")
	ntime <- 1
}
if (!exists("npls")){
	cat(" Number of PLS components not specified. \n Use argument 'npls= ' to specify the number of PLS components. \n Defaulting to 2 PLS components per trend.")
	npls <- 2
}
if (!exists("beta0cov")){
	cat(" Covariance form for beta0 not specified. \n Use argument 'beta0cov= ' to specify the covariance form. \n Defaulting to exponential.\n")
	beta0cov <- "exp"
}
if (!exists("beta_x_cov")){
	cat(" Covariance form for betaX (X = 1, 2, ..., ntime) not specified. \n Use argument 'beta_x_cov= ' to specify the covariance form. \n Defaulting to exponential.")
	beta_x_cov <- "exp"
}
if (!exists("nu_re")){
	cat(" Must specify whether to include random effect in nu field or not. Defaulting to including it.\n")
	nu_re <- TRUE
}
if (!exists("dftrend")){
	cat(" Must specify the number of degrees of freedom for the time trends. Defaulting to 152.\n")
	dftrend <- 152
}
if (! beta0cov %in% c("exp", "iid")){
	stop("Incorrect specification of the covariance form for beta0. \n Should be either 'exp' or 'iid'")
}
if (! beta_x_cov %in% c("exp", "iid")){
	stop("Incorrect specification of the covariance form for beta0. \n Should be either 'exp' or 'iid'")
}

#agencyFile <- paste("../dr0247/dr0247_L_agency_2005_2014_", pollutant, ".csv", sep = "")
#mesaFile <- paste("../dr0247/dr0247_L_mesa_PSD_summer_", pollutant, ".csv", sep = "")

##remove colocated with fixed
Agency<-Agency[which(Agency$native_id!="530530029"),]
Agency<-Agency[which(Agency$native_id!="530530031"),]
Agency<-Agency[which(Agency$native_id!="530332004"),]
Agency<-Agency[which(Agency$native_id!="530330057"),]
Agency<-Agency[which(Agency$native_id!="530330080"),]
Agency<-Agency[which(Agency$native_id!="530611007"),]
Agency<-Agency[which(Agency$native_id!="530610020"),]
dim(Agency)

##514   7
aqsData <- Agency

aqsData$PM25_concentration<-aqsData$concentration
aqsData$intended_wednesday<-aqsData$period_midpoint
aqsData$intended_wednesday<-as.Date(aqsData$intended_wednesday)


##
### Nancy will keep in colocated with fixed for CV 
##

##remove colocated with fixed
Low_cost<-Low_cost[which(Low_cost$native_id!="PSCA1"),]
Low_cost<-Low_cost[which(Low_cost$native_id!="PSCA2"),]
Low_cost<-Low_cost[which(Low_cost$native_id!="PSCA3"),]
Low_cost<-Low_cost[which(Low_cost$native_id!="PSCA4"),]
Low_cost<-Low_cost[which(Low_cost$native_id!="PSCA5"),]
Low_cost<-Low_cost[which(Low_cost$native_id!="PSCA6"),]

dim(Low_cost)
##[1] 184   7

Low_cost$site_type <- rep("C",length(Low_cost$concentration))
Fixed$site_type <- rep("F",length(Fixed$concentration))
mesaData <- rbind(Low_cost,Fixed)
mesaData$PM25_concentration<-mesaData$concentration
mesaData$intended_wednesday<-as.Date(mesaData$period_midpoint)
mesaData$mid<-mesaData$period_midpoint
mesaData$mid<-as.Date(mesaData$mid)
mesaData$day <- weekdays(as.Date(mesaData$mid))
mesaData$site_type <-as.factor(mesaData$site_type)

      #  mesaData$mid <- mesaData$start_day + floor((mesaData$end_day-mesaData$start_day)/2)
      #  mesaData$day <- weekdays(as.Date(mesaData$mid))
      #  mesaData<-mesaData[which(mesaData$day!="Saturday"),]
      #  dim(mesaData)
      #  mesaData<-mesaData[which(mesaData$day!="Sunday"),]
      #  dim(mesaData)
aqsConcColumn <- paste(pollutant, "_concentration", sep = "")
mesaConcColumn <- paste(pollutant, "_concentration", sep = "")

if (include_comco == FALSE) {
  mesaData <- mesaData[which(mesaData$site_type != "C"), ]
}



obs <- timeseries(aqsData, mesaData, aqsConcColumn, mesaConcColumn, npls)

covariatesRaw <- AQS_cov
covariatesRaw<-rbind(covariatesRaw,Low_cost_cov)
covariatesRaw<-rbind.fill(covariatesRaw,Fixed_cov)

covariatesRaw$region <- "Seattle"

preprocessedData <- covariatePrepare(covariatesRaw, obs, region) 

# can change the preset values based on the data you have 
range <- c(3.5, 4, 2, 8)
sill <- c( -3.2, -4, -2, -5)
nugget <- c(-3.2, -4, -2, -5)
random <- c(-3.5, -4, -2, -5)

#range <- c(10.33, 8, 5, 2.5)
#sill <- c( -4.09, -3, -2, -2)
#nugget <- c(-4.56, -3, -4.10, -7)
#random <- c(-2.607, -8, -5, -4)

modelResults <- model_fit_params(preprocessedData, npls=npls, ntime=ntime,
  beta0cov=beta0cov, beta_x_cov=beta_x_cov, nu_re=nu_re, cov_beta_nugget=cov_beta_nugget,
  cov_nu_nugget=cov_nu_nugget, range=range, sill=sill,
  nugget=nugget, random=random, dftrend=dftrend)

model.conv <- any(modelResults$st_estimate$summary$status$conv)

print("model convergence:")
print(modelResults$st_estimate$summary$status$conv)
print(modelResults$st_estimate$res.best)

plot(modelResults$st_data, "loc", main="Occurrence of Observations", xlab="", ylab="Location", col=c("black", "red","blue"), legend.loc=NULL)
#svg(filename="Occurrence_of_Observations.svg")
#plot(modelResults$st_data, "loc", main="Occurrence of Observations", xlab="", ylab="Location", col=c("black", "red"), legend.loc=NULL)
#dev.off()

# define all of the variables manually, try adding in the def for each of the variables within the function 

cvResults <- runCV(modelResults$st_model, 
                   modelResults$st_estimate, 
                   include_comco=include_comco, 
                   include_home=FALSE)

# describes the parameters for each cv 
# run again 


model_fname <- paste(outputDir, pollutant, ".", region, ".b0=", beta0cov,
 ".bx=", beta_x_cov, ".cnn=", cov_nu_nugget, ".cbn=", cov_beta_nugget, ".df=", dftrend,
 ".nure=", nu_re, ".nt=", ntime, ".npls=", npls,
 ".conv=", model.conv,".comco=",include_comco, ".Low_Cost.RData", sep = "")

print("saving to:")
print(model_fname)

# saves all of the model prediction results 
# CV results, prediction log scale, pred. LTA for long term, pred. obs for 2 week predictions 
# calculate R2 based on the observed and prediction 


save(list = c("obs", "covariatesRaw", "preprocessedData", "modelResults",
 "cvResults"), file = model_fname)

## additional function needed to summarize the results, for CV RMSE, R2, 


print("FINISHED")
