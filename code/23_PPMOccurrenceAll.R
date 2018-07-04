##This code estimates the percentage of posterior predictive mass greater or less than zero.
library(rstan)
library(bayesplot)
library(ggplot2)
library(gridExtra)
infolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/model_fits/"
outfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/summary_plots/"
load(paste0(infolder, "occ_fit_all_hillshade.RData")) #load model fit data
db.stan <- read.csv("/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/occurence_model_all_hillshade_noscale.csv")

#rename slope estimates (beta's)
names(occ_fit_all)[486:491] <- c("Elevation", "Hillshade", "Winter precip.",  "Winter precip./Total precip.", "Burned", "Grazed")

#Extract posterior samples for each slope parameter
posterior <- as.data.frame(occ_fit_all, pars=c("beta_p", "beta_o"))
#Percentage of posterior samples greater or less than zero
perc_pos <- apply(posterior, 2, function(x) (length(which(x > 0))/length(x) * 100))
perc_neg <- apply(posterior, 2, function(x) (length(which(x < 0))/length(x) * 100))

mass_sum <- cbind(perc_pos, perc_neg)
write.csv(mass_sum, file=paste0(outfolder,"Occ_All_Hillshade_MassSum.csv"))
