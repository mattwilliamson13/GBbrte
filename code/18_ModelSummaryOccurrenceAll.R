library(rstan)
library(bayesplot)
library(gridExtra)
infolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/model_fits/"
outfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/par_est_occ/"
load(paste0(infolder, "occ_fit_all_hillshade.RData")) #load model fit data
db.stan <- read.csv("/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/occurence_model_all_hillshade.csv")

#rename slope estimates (beta's)
names(occ_fit_all)[486:491] <- c("Elevation", "Hillshade", "Winter precip.", "Winter precip./Total precip.", "Burned", "Grazed")

#Create tabular slope results
s <- summary(occ_fit_all, pars=c("beta_p", "beta_o"), probs = c(0.10, 0.5, 0.90))
beta_sum <- as.data.frame(s$summary)
write.csv(beta_sum, file=paste0(outfolder,"occ_fit_all_hillshade_slopes.csv"))

#plot slope estimates
posterior <- as.array(occ_fit_all)
color_scheme_set("red")
slopes <- mcmc_intervals(posterior, regex_pars = c("Elevation", "Hillshade", "Winter precip.",  "Winter precip./Total precip.", "Burned", "Grazed"), prob=0.9)

slopes_area <- mcmc_areas(
  posterior,
  regex_pars =  c("Elevation", "Hillshade", "Winter precip.",  "Winter precip./Total precip.", "Burned", "Grazed"),
  prob = 0.9, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
)

g <- grid.arrange(slopes, slopes_area)
ggsave(g, filename = paste0(outfolder,"occ_fit_all_hillshade_slopes.tiff"))

##INTERCEPTS
#rename range intercepts
range_names <- levels(db.stan$range)
names(occ_fit_all)[492:495] <- paste0("a_",range_names)
#rename canyon intercepts
can_names <- levels(db.stan$area)
names(occ_fit_all)[496:524] <- paste0("a_",can_names)

#create lookup table for Can in Rg
cr <- db.stan[,c(1:2, 13:14)]
cr.u <- unique(cr[,c('range','area','area_id','range_id')])

cans.in.r.1 <- paste0("a_",as.character(cr.u[cr.u$range_id==1, 2]))
cans.in.r.2 <- paste0("a_",as.character(cr.u[cr.u$range_id==2, 2]))
cans.in.r.3 <- paste0("a_",as.character(cr.u[cr.u$range_id==3, 2]))
cans.in.r.4 <- paste0("a_",as.character(cr.u[cr.u$range_id==4, 2]))

#Create tabular range intercept results
s <- summary(occ_fit_all, pars=c("mu_alpha","alpha_rg"), probs = c(0.10, 0.5, 0.90))
alpha_rg_sum <- as.data.frame(s$summary)
write.csv(alpha_rg_sum, file=paste0(outfolder,"occ_fit_all_hillshade_int_rg.csv"))

#Create tabular canyon intercept results
s <- summary(occ_fit_all, pars=c("alpha_can"), probs = c(0.10, 0.5, 0.90))
alpha_can_sum <- as.data.frame(s$summary)
write.csv(alpha_can_sum, file=paste0(outfolder,"occ_fit_all_hillshade_int_can.csv"))



#plot range intercept estimates
posterior <- as.array(occ_fit_all)
color_scheme_set("red")
int_rg.1 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[1]), cans.in.r.1), prob=0.90) + ggplot2::labs(title=range_names[1])
int_rg.2 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[2]), cans.in.r.2), prob=0.90) + ggplot2::labs(title=range_names[2])
int_rg.3 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[3]), cans.in.r.3), prob=0.90) + ggplot2::labs(title=range_names[3])
int_rg.4 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[4]), cans.in.r.4), prob=0.90) + ggplot2::labs(title=range_names[4])

int_rg <- grid.arrange(int_rg.1, int_rg.2, int_rg.3, int_rg.4)
ggsave(int_rg, filename = paste0(outfolder,"occ_fit_all_hillshade_ints.tiff"))

int_rg_areas.1 <- mcmc_areas(
  posterior,
  regex_pars =  c("mu_alpha",paste0("a_", range_names[1]), cans.in.r.1),
  prob = 0.90, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[1])
int_rg_areas.2 <- mcmc_areas(
  posterior,
  regex_pars = c("mu_alpha",paste0("a_", range_names[2]), cans.in.r.2),
  prob = 0.90, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[2])
int_rg_areas.3 <- mcmc_areas(
  posterior,
  regex_pars =  c("mu_alpha",paste0("a_", range_names[3]), cans.in.r.3),
  prob = 0.90, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[3])
int_rg_areas.4 <- mcmc_areas(
  posterior,
  regex_pars =  c("mu_alpha",paste0("a_", range_names[4]), cans.in.r.4),
  prob = 0.90, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[4])

int_rg_ar<- grid.arrange(int_rg_areas.1, int_rg_areas.2, int_rg_areas.3, int_rg_areas.4)
ggsave(int_rg_ar, filename = paste0(outfolder,"occ_fit_all_hillshade_ints_area.tiff"))
