#estimate the posterior predictive mass for regression parameters in the occurrence model for unburned data
library(rstan)
library(bayesplot)
library(ggplot2)
library(gridExtra)
infolder <- here::here("model_fits/")
outfolder <- here::here("param_est/")
load(paste0(infolder, "/occ_fit_ub_hillshade.RData")) #load model fit data
db.stan <- read.csv(here::here("processed_data/occurence_model_ub_hillshade_noscale.csv"))

#rename slope estimates (beta's)
names(occ_fit_ub)[391:395] <- c("Elevation",  "Hillshade", "Winter precip.", "Winter precip./Total precip.", "Grazed")

#Extract posterior samples for each slope parameter
posterior <- as.data.frame(occ_fit_ub, pars=c("beta_p", "beta_o"))
#Percentage of posterior samples greater or less than zero
perc_pos <- apply(posterior, 2, function(x) (length(which(x > 0))/length(x) * 100))
perc_neg <- apply(posterior, 2, function(x) (length(which(x < 0))/length(x) * 100))

mass_sum <- cbind(perc_pos, perc_neg)
write.csv(mass_sum, file=paste0(outfolder,"/Occ_UB_hillshade_MassSum.csv"))
