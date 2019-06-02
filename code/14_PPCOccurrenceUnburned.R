library(pROC)
library(rstan)
library(LaplacesDemon)
infolder <- here::here("model_fits/")
load(paste0(infolder, "/occ_fit_ub_hillshade.RData")) #load model fit data
db.stan <- read.csv(here::here("processed_data/occurence_model_unburned_hillshade.csv"))
thetas <- summary(occ_fit_ub, pars='theta', probs = c(0.25,0.5,0.75)) #summarize posterior estimates of theta

#extract summaries of posterior probability of occurrence
th <- as.data.frame(thetas$summary)
th.med <- th$`50%` #get median value of posterior theta estimate
th.25 <- th$`25%` #get lower quartile value of posterior theta estimate
th.75 <- th$`75%` #get upper quartile value of posterior theta estimate

#generate rocs based on summary stats of posterior
roc(db.stan$occ, invlogit(th.25), plot = TRUE, main = "ROC for occurrence data on unburned sites \n based on the lower quartile estimates of theta") #AUC = .912
roc(db.stan$occ, invlogit(th.med), plot=TRUE, main = "ROC for occurrence data on unburned sites \n based on the median estimates of theta") #AUC = .9135
roc(db.stan$occ, invlogit(th.75), plot = TRUE, main = "ROC for occurrence data on unburned sites \n based on the upper quartile estimates of theta") # AUC =.914
