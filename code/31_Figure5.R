library(extrafont)
library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(rethinking)
library(egg)
library(viridis)
library(tidyverse)
infolder <- here::here("model_fits/")
outfolder <- here::here("summary_plots/")
datafolder <- here::here("processed_data/")

#load the two occurrence data sets
occ.all.ns <- read.csv(paste0(datafolder,"/occurence_model_all_hillshade_noscale.csv"), stringsAsFactors = FALSE)

#boxplot of the actual distribution of BRTE occurrence across elevation
elev.comb <- ggplot(occ.all.ns, aes(x = "All points", y=elev, fill=as.factor(occ))) +
  geom_boxplot() +
  geom_boxplot(data=occ.all.ns, aes(x=as.factor(burn), y=elev, fill=as.factor(occ))) +
  scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL) +
  scale_x_discrete(name = NULL, labels = c("All \npoints", "Burned \npoints", "Unburned \npoints")) +
  labs(y="Elevation (m)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

#Marginal effect of elevation effect on occurrence
load(paste0(infolder, "/occ_fit_all_hillshade.RData"))
posterior.occ <- as.data.frame(occ_fit_all, pars=c("beta_p", "beta_o"))
#get the intercepts
int.occ <- as.matrix(occ_fit_all, pars="mu_alpha")
#extract the burn covariate
post_b.occ <- as.matrix(posterior.occ[,5])
#extract the elevation effect
post_e.occ <- as.matrix(posterior.occ[,1])

occ.all <- read.csv(paste0(datafolder,"/occurence_model_all_hillshade.csv"), stringsAsFactors = FALSE)

elev.occ.b = as.matrix(seq(min(occ.all[occ.all$burn== 1,]$elev), max(occ.all[occ.all$burn== 1,]$elev), length.out=nrow(occ.all[occ.all$burn== 1,])))
elev.occ.ub = as.matrix(seq(min(occ.all[occ.all$burn== 0,]$elev), max(occ.all[occ.all$burn== 0,]$elev), length.out=nrow(occ.all[occ.all$burn== 0,])))
b.occ = as.matrix(rep(1,nrow(elev.occ.b)))
ub.occ = as.matrix(rep(0,nrow(elev.occ.ub)))

#estimate the effect
bep.b <- (post_e.occ %*% t(elev.occ.b)) + (post_b.occ %*% t(b.occ))
bep.ub <- (post_e.occ %*% t(elev.occ.ub)) + (post_b.occ %*% t(ub.occ))

marg.elev.b <- logistic(as.vector(int.occ) + bep.b)
marg.elev.ub <- logistic(as.vector(int.occ) + bep.ub)

probs_elev.b <- melt(marg.elev.b, value.name = "prob", varnames = c("iter", "elevid"))
probs_elev.b$burn <- 1

probs_elev.ub <- melt(marg.elev.ub, value.name = "prob", varnames = c("iter", "elevid"))
probs_elev.ub$burn <- 0

elev.b_df <- data.frame(
  elev = as.vector(elev.occ.b),
  id = seq(1:length(elev.occ.b))
)

elev.ub_df <- data.frame(
  elev = as.vector(elev.occ.ub),
  id = seq(1:length(elev.occ.ub))
)

probs_elev.b$elev <- elev.b_df$elev[probs_elev.b$elevid]
probs_elev.b$elevconv <- (probs_elev.b$elev * sd(occ.all.ns$elev)) + mean(occ.all.ns$elev)

probs_elev.ub$elev <- elev.ub_df$elev[probs_elev.ub$elevid]
probs_elev.ub$elevconv <- (probs_elev.ub$elev * sd(occ.all.ns$elev)) + mean(occ.all.ns$elev)

probs_elev <- rbind(probs_elev.b, probs_elev.ub)

prob_elev.occ <- probs_elev %>%
  group_by(elevconv,burn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )

p_elev.occ <- ggplot(prob_elev.occ) + geom_line(aes(y=p, x=elevconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=elevconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Elevation (m)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))


#Marginal effects of elevation on prevalence for all three models (fit to all data, unburned only data, and burned only data)
#load the prevalence model fits
load(paste0(infolder, "/reprate_fit_all_hillshade.RData"))
load(paste0(infolder, "/reprate_fit_ub_hillshade.RData")) #load model fit data
load(paste0(infolder, "/reprate_fit_b_hillshade2.RData")) #load model fit data


posterior.rr.all <- as.data.frame(reprate_fit_all, pars=c("beta_p", "beta_o"))
posterior.rr.b <- as.data.frame(reprate_fit_b, pars=c("beta_p", "beta_o"))
posterior.rr.ub <- as.data.frame(reprate_fit_ub, pars=c("beta_p", "beta_o"))

#get the intercepts
int.rr.all <- as.matrix(reprate_fit_all, pars="mu_alpha")
int.rr.b <- as.matrix(reprate_fit_b, pars="mu_alpha")
int.rr.ub <- as.matrix(reprate_fit_ub, pars="mu_alpha")

#extract the burn covariate
post_b.rr.all <- as.matrix(posterior.rr.all[,5])

#extract the elevation covariate
post_e.rr.all <- as.matrix(posterior.rr.all[,1])
post_e.rr.ub <- as.matrix(posterior.rr.ub[,1])
post_e.rr.b <- as.matrix(posterior.rr.b[,1])

#set up estimation matrices
rr.all <- read.csv(paste0(datafolder,"/reprate_model_all_hillshade.csv"), stringsAsFactors = FALSE)
rr.b <- read.csv(paste0(datafolder,"/reprate_model_burned_hillshade.csv"), stringsAsFactors = FALSE)
rr.ub <- read.csv(paste0(datafolder,"/reprate_model_unburned_hillshade.csv"), stringsAsFactors = FALSE)

elev.rr.b.all = as.matrix(seq(min(rr.all[rr.all$burn== 1,]$elev), max(rr.all[rr.all$burn== 1,]$elev), length.out=nrow(rr.all[rr.all$burn== 1,])))
elev.rr.ub.all = as.matrix(seq(min(rr.all[rr.all$burn== 0,]$elev), max(rr.all[rr.all$burn== 0,]$elev), length.out=nrow(rr.all[rr.all$burn== 0,])))
elev.rr.b <- as.matrix(seq(min(rr.b$elev), max(rr.b$elev), length.out=nrow(rr.b)))
elev.rr.ub <- as.matrix(seq(min(rr.ub$elev), max(rr.ub$elev), length.out=nrow(rr.ub)))
b.rr = as.matrix(rep(1,nrow(elev.rr.b.all)))
ub.rr = as.matrix(rep(0,nrow(elev.rr.ub.all)))

#estimate effect

bwe.b.all <- (post_e.rr.all %*% t(elev.rr.b.all)) + (post_b.rr.all %*% t(b.rr))
bwe.ub.all <- (post_e.rr.all %*% t(elev.rr.ub.all)) + (post_b.rr.all %*% t(ub.rr))
bwe.b <- (post_e.rr.b %*% t(elev.rr.b))
bwe.ub <- (post_e.rr.ub %*% t(elev.rr.ub))


#generate posterior predictions
marg.elev.b.all <- logistic(as.vector(int.rr.all) + bwe.b.all)
marg.elev.ub.all <- logistic(as.vector(int.rr.all) + bwe.ub.all)
marg.elev.b <- logistic(as.vector(int.rr.b) + bwe.b)
marg.elev.ub <- logistic(as.vector(int.rr.ub) + bwe.ub)

#create dataframe
probs_elev_b.all <- melt(marg.elev.b.all, value.name = "prob", varnames = c("iter", "elevid"))
probs_elev_b.all$burn <- 1
probs_elev_ub.all <- melt(marg.elev.ub.all, value.name = "prob", varnames = c("iter", "elevid"))
probs_elev_ub.all$burn <- 0

probs_elev_b <- melt(marg.elev.b, value.name = "prob", varnames = c("iter", "elevid"))
probs_elev_b$burn <- 1
probs_elev_ub <- melt(marg.elev.ub, value.name = "prob", varnames = c("iter", "elevid"))
probs_elev_ub$burn <- 0

elev.b.all_df <- data.frame(
  elev = as.vector(elev.rr.b.all),
  id = seq(1:length(elev.rr.b.all))
)

elev.ub.all_df <- data.frame(
  elev = as.vector(elev.rr.ub.all),
  id = seq(1:length(elev.rr.ub.all))
)


elev.b_df <- data.frame(
  elev = as.vector(elev.rr.b),
  id = seq(1:length(elev.rr.b))
)

elev.ub_df <- data.frame(
  elev = as.vector(elev.rr.ub),
  id = seq(1:length(elev.rr.ub))
)

# join grazing data to probs_df
rr.all.ns <- read.csv(here::here("/processed_data/reprate_model_all_noscale_hillshade.csv"))
probs_elev_b.all$elev <- elev.b.all_df$elev[probs_elev_b.all$elevid]
probs_elev_b.all$elevconv <- (probs_elev_b.all$elev * sd(rr.all.ns$elev)) + mean(rr.all.ns$elev)
probs_elev_ub.all$elev <- elev.ub.all_df$elev[probs_elev_ub.all$elevid]
probs_elev_ub.all$elevconv <- (probs_elev_ub.all$elev * sd(rr.all.ns$elev)) + mean(rr.all.ns$elev)

probs_df.all <- rbind(probs_elev_b.all, probs_elev_ub.all)

rr.b.ns <- read.csv(here::here("/processed_data/reprate_model_burned_noscale_hillshade.csv"))
probs_elev_b$elev <- elev.b_df$elev[probs_elev_b$elevid]
probs_elev_b$elevconv <- (probs_elev_b$elev * sd(rr.b.ns$elev)) + mean(rr.b.ns$elev)

rr.ub.ns <- read.csv(here::here("/processed_data/reprate_model_unburned_noscale_hillshade.csv"))
probs_elev_ub$elev <- elev.ub_df$elev[probs_elev_ub$elevid]
probs_elev_ub$elevconv <- (probs_elev_ub$elev * sd(rr.ub.ns$elev)) + mean(rr.ub.ns$elev)

probs_df.comb <- rbind(probs_elev_b, probs_elev_ub)



# compute summaries and then plot
prob_e.all <- probs_df.all %>%
  group_by(burn,elevconv) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )

prob_e.comb <- probs_df.comb %>%
  group_by(burn,elevconv) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )


p_elev.all <- ggplot(prob_e.all) + geom_line(aes(y=p, x=elevconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=elevconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Elevation (m)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

p_elev.comb <- ggplot(prob_e.comb) + geom_line(aes(y=p, x=elevconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=elevconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Elevation (SD)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))


elev.plot <- plot_grid(elev.comb,
                       p_elev.occ + ylim(c(0,1)),
                       p_elev.all + ylim(c(0,1)),
                       p_elev.comb + ylab("") + ylim(c(0,1)),
                       align="hv",
                       axis="bl",
                       nrow=2,
                       labels = c("A", "B","C","D"),
                       label_x = c(0.89, 0.89, 0.89, 0.89),
                       label_y = c(0.98,0.98,0.98, 0.98),
                       label_fontfamily = "Times New Roman",
                       vjust = 1
                      )
cowplot::ggsave(paste0(outfolder,"/Figure5.tiff"),plot=elev.plot, height= 6, width=6, units="in")
