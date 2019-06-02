library(pROC)
library(rstan)
library(LaplacesDemon)
infolder <- here::here("model_fits/")
load(paste0(infolder, "/occ_fit_all_hillshade.RData")) #load model fit data
db.stan <- read.csv(here::here("/processed_data/occurence_model_all_hillshade.csv"))
thetas <- summary(occ_fit_all, pars='theta', probs = c(0.25,0.5,0.75)) #summarize posterior estimates of theta

#extract summary statistics of posterior probability
th <- as.data.frame(thetas$summary)
th.med <- th$`50%` #get median value of posterior theta estimate
th.25 <- th$`25%` #get lower quartile value of posterior theta estimate
th.75 <- th$`75%` #get upper quartile value of posterior theta estimate

#fit ROC based on posterior summaries
roc(db.stan$occ, invlogit(th.25), plot = TRUE, main = "ROC for all occurrence data based on the lower quartile estimates of theta") #AUC = .935
roc(db.stan$occ, invlogit(th.med), plot=TRUE, main = "ROC for all occurrence data based on the median estimates of theta") #AUC = .936
roc(db.stan$occ, invlogit(th.75), plot = TRUE, main = "ROC for all occurrence data based on the upper quartile estimates of theta") # AUC =.9368
