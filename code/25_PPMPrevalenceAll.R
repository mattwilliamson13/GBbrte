#estimate the posterior predictive mass for predictors in the prevalence models for all points

library(rstan)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(egg)
library(rethinking)
infolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/model_fits/"
outfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/summary_plots/"
load(paste0(infolder, "reprate_fit_all_hillshade.RData")) #load model fit data
db.stan <- read.csv("/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/reprate_model_all_hillshade.csv")

#rename slope estimates (beta's)
names(reprate_fit_all)[334:343] <- c("Elevation","Hillshade","Med. Winter Precip.", "Med. Winter precip./Total precip.", "Burned", "Grazed", "Perennial grass frequency", "Spring precip.", "Winter precip.", "Winter precip./Total precip.")

#Extract posterior samples for each slope parameter
posterior <- as.data.frame(reprate_fit_all, pars=c("beta_p", "beta_o"))
#Percentage of posterior samples greater or less than zero
perc_pos <- apply(posterior, 2, function(x) (length(which(x > 0))/length(x) * 100))
perc_neg <- apply(posterior, 2, function(x) (length(which(x < 0))/length(x) * 100))

mass_sum <- cbind(perc_pos, perc_neg)
write.csv(mass_sum, file=paste0(outfolder,"RR_All_MassSum_hillshade.csv"))
