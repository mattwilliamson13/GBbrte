library(loo)
library(bayesplot)
options(mc.cores = 6)
infolder <- here::here("model_fits/")

load(paste0(infolder, "/reprate_fit_all_hillshade.RData")) #load model fit data
db.stan <- read.csv(here::here("processed_data/reprate_model_all_hillshade.csv"), stringsAsFactors = FALSE)

###density plot overlay
y_ppcs <- as.matrix(reprate_fit_all, pars = "y_tilde")
ppc_dens_overlay(
  y = db.stan$brte,
  yrep = y_ppcs[sample(nrow(y_ppcs), 100),]
)

###Test statistic graphics
ppc_stat(
  y = db.stan$brte,
  yrep = y_ppcs,
  stat="mean",
  binwidth = 0.05
)


ppc_stat(
  y = db.stan$brte,
  yrep = y_ppcs,
  stat="median",
  binwidth = 0.05
)

ppc_stat(
  y = db.stan$brte,
  yrep = y_ppcs,
  stat="sd",
  binwidth = 0.05
)

ppc_stat(
  y = db.stan$brte,
  yrep = y_ppcs,
  stat="max",
  binwidth = 0.05
)

##grouped by burn
color_scheme_set("teal")
ppc_stat_grouped(
  y = db.stan$brte,
  yrep = y_ppcs,
  group = db.stan$burn,
  stat="mean",
  binwidth = 0.05
)

ppc_stat_grouped(
  y = db.stan$brte,
  yrep = y_ppcs,
  group = db.stan$burn,
  stat="median",
  binwidth = 0.05
)

ppc_stat_grouped(
  y = db.stan$brte,
  yrep = y_ppcs,
  group = db.stan$burn,
  stat="sd",
  binwidth = 0.05
)

ppc_stat_grouped(
  y = db.stan$brte,
  yrep = y_ppcs,
  group = db.stan$burn,
  stat="max",
  binwidth = 0.05
)

##Numerical fit estimates
rr_sim <- rstan::extract(reprate_fit_all)

od2 <- rr_sim$obs_dev2
rd2 <- rr_sim$sim_dev2
fit <- apply(od2, 2, sum)
fit.new <- apply(rd2, 2, sum)
gof <- mean(ifelse((fit.new - fit) > 0, 1, 0)) #0.2564

##LOO estimate provides an estimate of the number of points that are extreme in light of the model
log_lik_fit <-  extract_log_lik(reprate_fit_all, parameter_name = "log_lik")

loo_fit <- loo(log_lik_fit)

