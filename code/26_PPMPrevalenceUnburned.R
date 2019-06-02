#estimate the posterior predictive mass for predictors in models of prevalence applied to the unburned points

library(rstan)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(egg)
library(rethinking)
infolder <- here::here("model_fits/")
outfolder <- here::here("param_est/")
load(paste0(infolder, "/reprate_fit_ub_hillshade.RData")) #load model fit data
db.stan <- read.csv(here::here("processed_data/reprate_model_unburned_hillshade.csv"))

#rename slope estimates (beta's)
names(reprate_fit_ub)[247:254] <- c("Elevation", "Hillshade","Med. Winter precip.", "Med. Winter precip./Total precip." ,"Grazed", "Perennial grass frequency", "Spring precip.", "Winter precip./Total precip.")

#Extract posterior samples for each slope parameter
posterior <- as.data.frame(reprate_fit_ub, pars=c("beta_p", "beta_o"))
#Percentage of posterior samples greater or less than zero
perc_pos <- apply(posterior, 2, function(x) (length(which(x > 0))/length(x) * 100))
perc_neg <- apply(posterior, 2, function(x) (length(which(x < 0))/length(x) * 100))

mass_sum <- cbind(perc_pos, perc_neg)
write.csv(mass_sum, file=paste0(outfolder,"/RR_UB_MassSum_hillshade.csv"))
