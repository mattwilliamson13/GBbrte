library(rstan)
library(bayesplot)
library(gridExtra)
library(gridExtra)
infolder <- here::here("model_fits/")
outfolder <- here::here("param_est/")
load(paste0(infolder, "/reprate_fit_ub_hillshade.RData")) #load model fit data
db.stan <- read.csv(here::here("processed_data/reprate_model_unburned_hillshade.csv"))

#rename slope estimates (beta's)
names(reprate_fit_ub)[247:254] <- c("Elevation", "Hillshade","Med. Winter precip.", "Med. Winter precip./Total precip." ,"Grazed", "Perennial grass frequency", "Spring precip.", "Winter precip./Total precip.")

#Create tabular slope results
s <- summary(reprate_fit_ub, pars=c("beta_p", "beta_o"), probs = c(0.1, 0.5, 0.9))
beta_sum <- as.data.frame(s$summary)
write.csv(beta_sum, file=paste0(outfolder,"/reprate_fit_ub_slopes_hillshade.csv"))

#plot slope estimates
posterior <- as.array(reprate_fit_ub)
color_scheme_set("red")
slopes <- mcmc_intervals(posterior, regex_pars = c("Elevation", "Hillshade","Med. Winter precip.", "Med. Winter precip./Total precip.", "Grazed", "Perennial grass frequency", "Spring precip.", "Winter precip./Total precip."), prob=0.90)

slopes_area <- mcmc_areas(
  posterior,
  regex_pars =  c("Elevation", "Hillshade","Med. Winter precip.", "Med. Winter precip./Total precip.", "Grazed", "Perennial grass frequency", "Spring precip.", "Winter precip./Total precip."),
  prob = 0.9, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
)

g <- grid.arrange(slopes, slopes_area)
ggsave(g, filename = paste0(outfolder,"/reprate_fit_ub_slopes_hillshade.tiff"))

#INTERCEPTS
#rename range intercepts
range_names <- levels(db.stan$range)
names(reprate_fit_ub)[255:258] <- paste0("a_",range_names)
#rename canyon intercepts
can_names <- levels(db.stan$area)
names(reprate_fit_ub)[259:283] <- paste0("a_",can_names)

#create lookup table for Can in Rg
cr <- db.stan[,c(1:2, 14:15)]
cr.u <- unique(cr[,c('range','area','area_id','range_id')])

cans.in.r.1 <- paste0("a_",as.character(cr.u[cr.u$range_id==1, 2]))
cans.in.r.2 <- paste0("a_",as.character(cr.u[cr.u$range_id==2, 2]))
cans.in.r.3 <- paste0("a_",as.character(cr.u[cr.u$range_id==3, 2]))
cans.in.r.4 <- paste0("a_",as.character(cr.u[cr.u$range_id==4, 2]))

#Create tabular range intercept results
s <- summary(reprate_fit_ub, pars=c("alpha_rg"), probs = c(0.1, 0.5, 0.9))
alpha_rg_sum <- as.data.frame(s$summary)
write.csv(alpha_rg_sum, file=paste0(outfolder,"/reprate_fit_ub_int_rg_hillshade.csv"))

#Create tabular canyon intercept results
s <- summary(reprate_fit_ub, pars=c("alpha_can"), probs = c(0.1, 0.5, 0.9))
alpha_can_sum <- as.data.frame(s$summary)
write.csv(alpha_can_sum, file=paste0(outfolder,"/reprate_fit_ub_int_can_hillshade.csv"))

#plot range intercept estimates
posterior <- as.array(reprate_fit_ub)
color_scheme_set("red")
int_rg.1 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[1]), cans.in.r.1), prob=0.9) + ggplot2::labs(title=range_names[1])
int_rg.2 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[2]), cans.in.r.2), prob=0.9) + ggplot2::labs(title=range_names[2])
int_rg.3 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[3]), cans.in.r.3), prob=0.9) + ggplot2::labs(title=range_names[3])
int_rg.4 <- mcmc_intervals(posterior, regex_pars = c("mu_alpha",paste0("a_", range_names[4]), cans.in.r.4), prob=0.9) + ggplot2::labs(title=range_names[4])

int_rg <- grid.arrange(int_rg.1, int_rg.2, int_rg.3, int_rg.4)
ggsave(int_rg, filename = paste0(outfolder,"/reprate_fit_ub_ints_hillshade.tiff"))

int_rg_areas.1 <- mcmc_areas(
  posterior,
  regex_pars =  c("mu_alpha",paste0("a_", range_names[1]), cans.in.r.1),
  prob = 0.9, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[1])
int_rg_areas.2 <- mcmc_areas(
  posterior,
  regex_pars = c("mu_alpha",paste0("a_", range_names[2]), cans.in.r.2),
  prob = 0.9, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[2])
int_rg_areas.3 <- mcmc_areas(
  posterior,
  regex_pars =  c("mu_alpha",paste0("a_", range_names[3]), cans.in.r.3),
  prob = 0.9, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[3])
int_rg_areas.4 <- mcmc_areas(
  posterior,
  regex_pars =  c("mu_alpha",paste0("a_", range_names[3]), cans.in.r.4),
  prob = 0.9, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) + ggplot2::labs(title=range_names[4])

int_rg_ar<- grid.arrange(int_rg_areas.1, int_rg_areas.2, int_rg_areas.3, int_rg_areas.4)
ggsave(int_rg_ar, filename = paste0(outfolder,"/reprate_fit_ub_ints_area_hillshade.tiff"))
