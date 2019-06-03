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
occ.all <- read.csv(paste0(datafolder,"/occurence_model_all_hillshade_noscale.csv"), stringsAsFactors = FALSE)
#generate boxplots of observed data for various values of winter precip
wpcp.comb <- ggplot(occ.all, aes(x = "All points", y=medwtpcp, fill=as.factor(occ))) +
  geom_boxplot() +
  geom_boxplot(data=occ.all, aes(x=as.factor(burn), y=medwtpcp, fill=as.factor(occ))) +
  scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL) +
  scale_x_discrete(name = NULL, labels = c("All \npoints", "Burned \npoints", "Unburned \npoints")) +
  labs(y="Median winter \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))
#generate boxplots of observed data for various values of spring precip
spcp.comb <- ggplot(occ.all, aes(x = "All points", y=medsppcp, fill=as.factor(occ))) +
  geom_boxplot() +
  geom_boxplot(data=occ.all, aes(x=as.factor(burn), y=medwtpcp, fill=as.factor(occ))) +
  scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL) +
  scale_x_discrete(name = NULL, labels = c("All \npoints", "Burned \npoints", "Unburned \npoints")) +
  labs(y="Median spring \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

#generate boxplots of observed data for various values of proportion of precip falling in winter
propwtr.comb <- ggplot(occ.all, aes(x = "All points", y=medpropwtr, fill=as.factor(occ))) +
  geom_boxplot() +
  geom_boxplot(data=occ.all, aes(x=as.factor(burn), y=medpropwtr, fill=as.factor(occ))) +
  scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL) +
  scale_x_discrete(name = NULL, labels = c("All \npoints", "Burned \npoints", "Unburned \npoints")) +
  labs(y="Median winter precipitation:\ntotal precipitation") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

pcp.plot.1 <- plot_grid(wpcp.comb,
                        spcp.comb,
                        propwtr.comb,
                        nrow=1,
                        align = "h",
                        axis="b")

##marginal effects of winter precipitation
#load the model fits
load(paste0(infolder, "/reprate_fit_all_hillshade.RData"))
load(paste0(infolder, "/occ_fit_all_hillshade.RData"))

posterior.rr <- as.data.frame(reprate_fit_all, pars=c("beta_p", "beta_o"))
posterior.occ <- as.data.frame(occ_fit_all, pars=c("beta_p", "beta_o"))

#get the intercepts
int.rr <- as.matrix(reprate_fit_all, pars="mu_alpha")
int.occ <- as.matrix(occ_fit_all, pars="mu_alpha")

#extract the burn covariate
post_b.rr <- as.matrix(posterior.rr[,5])
post_b.occ <- as.matrix(posterior.occ[,5])



######Precip effect on occurrence
#load scaled values
occ.all <- read.csv(paste0(datafolder,"/occurence_model_all_hillshade.csv"), stringsAsFactors = FALSE)

#extract winter precip cov
post_wp.occ <- as.matrix(posterior.occ[,3])
#extract proportion cov
post_pwp.occ <- as.matrix(posterior.occ[,4])

##NOTE: Have to use actual spring precip values here to determine appropriate values for the proportion
wtpcp.occ.b = as.matrix(seq(min(occ.all[occ.all$burn== 1,]$medwtpcp), max(occ.all[occ.all$burn== 1,]$medwtpcp), length.out=nrow(occ.all[occ.all$burn== 1,])))
wtpcp.occ.ub = as.matrix(seq(min(occ.all[occ.all$burn== 0,]$medwtpcp), max(occ.all[occ.all$burn== 0,]$medwtpcp), length.out=nrow(occ.all[occ.all$burn== 0,])))
b.occ = as.matrix(rep(1,nrow(wtpcp.occ.b)))
ub.occ = as.matrix(rep(0,nrow(wtpcp.occ.ub)))
sppcp.occ.b = as.matrix(rep(max(occ.all[occ.all$burn== 1,]$medsppcp), length.out=nrow(occ.all[occ.all$burn== 1,])))
sppcp.occ.ub = as.matrix(rep(max(occ.all[occ.all$burn== 0,]$medsppcp), length.out=nrow(occ.all[occ.all$burn== 0,])))


#load unscaled data to get proportion correct
occ.us <- read.csv(paste0(datafolder,"/occurence_model_all_hillshade_noscale.csv"), stringsAsFactors = FALSE)
#generate transformed values
wtpcp.occ.b.us <- wtpcp.occ.b * sd(occ.us$medwtpcp) + mean(occ.us$medwtpcp) #convert back to original scale
wtpcp.occ.ub.us <- wtpcp.occ.ub * sd(occ.us$medwtpcp) + mean(occ.us$medwtpcp) #convert back to original scale
sppcp.occ.b.us <- sppcp.occ.b * sd(occ.us$medsppcp) + mean(occ.us$medsppcp)
sppcp.occ.ub.us <- sppcp.occ.ub * sd(occ.us$medsppcp) + mean(occ.us$medsppcp)
#generate proportion and rescale
propwtr.b = scale(as.matrix(wtpcp.occ.b.us/(wtpcp.occ.b.us + sppcp.occ.b.us)))
propwtr.ub = scale(as.matrix(wtpcp.occ.ub.us/(wtpcp.occ.ub.us + sppcp.occ.ub.us)))


bwp.b <- (post_wp.occ %*% t(wtpcp.occ.b)) + (post_b.occ %*% t(b.occ)) +
          (post_pwp.occ %*% t(propwtr.b))
bwp.ub <- (post_wp.occ %*% t(wtpcp.occ.ub)) + (post_b.occ %*% t(ub.occ)) +
  (post_pwp.occ %*% t(propwtr.ub))


marg.wpcp.b.h <- logistic(as.vector(int.occ) + bwp.b)
marg.wpcp.ub.h <- logistic(as.vector(int.occ) + bwp.ub)

probs_wpcp.b.h <- melt(marg.wpcp.b.h, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.b.h$burn <- 1

probs_wpcp.ub.h <- melt(marg.wpcp.ub.h, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.ub.h$burn <- 0

wpcp.b_df <- data.frame(
  wpcp = as.vector(wtpcp.occ.b),
  id = seq(1:length(wtpcp.occ.b))
)

wpcp.ub_df <- data.frame(
  wpcp = as.vector(wtpcp.occ.ub),
  id = seq(1:length(wtpcp.occ.ub))
)

probs_wpcp.b.h$wpcp <- wpcp.b_df$wpcp[probs_wpcp.b.h$wpid]
probs_wpcp.b.h$wpcpconv <- (probs_wpcp.b.h$wpcp * sd(occ.us$medwtpcp)) + mean(occ.us$medwtpcp)

probs_wpcp.ub.h$wpcp <- wpcp.ub_df$wpcp[probs_wpcp.ub.h$wpid]
probs_wpcp.ub.h$wpcpconv <- (probs_wpcp.ub.h$wpcp * sd(occ.us$medwtpcp)) + mean(occ.us$medwtpcp)

probs_wpcp.h <- rbind(probs_wpcp.b.h, probs_wpcp.ub.h)

prob_wp.occ.high <- probs_wpcp.h %>%
  group_by(wpcpconv,burn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )

##Marginal effects of winter precip on occurrence (while holding spring precip at its maximum)
p_wpcp.occ.h <- ggplot(prob_wp.occ.high) + geom_line(aes(y=p, x=wpcpconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=wpcpconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Median winter \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

##Marginal effects of winter precip on occurrence (while holding spring precip at its mean)
#generate values for the mean spring precip
sppcp.occ.b.mean = as.matrix(rep(mean(occ.all[occ.all$burn== 1,]$medsppcp), length.out=nrow(occ.all[occ.all$burn== 1,])))
sppcp.occ.ub.mean = as.matrix(rep(mean(occ.all[occ.all$burn== 0,]$medsppcp), length.out=nrow(occ.all[occ.all$burn== 0,])))

#transform back to original scales for proportion
sppcp.occ.b.mean.us <- sppcp.occ.b.mean * sd(occ.us$medsppcp) + mean(occ.us$medsppcp)
sppcp.occ.ub.mean.us <- sppcp.occ.ub.mean * sd(occ.us$medsppcp) + mean(occ.us$medsppcp)

#generate proportion and rescale
propwtr.b.mean = scale(as.matrix(wtpcp.occ.b.us/(wtpcp.occ.b.us + sppcp.occ.b.mean.us)))
propwtr.ub.mean = scale(as.matrix(wtpcp.occ.ub.us/(wtpcp.occ.ub.us + sppcp.occ.ub.mean.us)))

bwp.b.mean <- (post_wp.occ %*% t(wtpcp.occ.b)) + (post_b.occ %*% t(b.occ)) +
  (post_pwp.occ %*% t(propwtr.b.mean))
bwp.ub.mean <- (post_wp.occ %*% t(wtpcp.occ.ub)) + (post_b.occ %*% t(ub.occ)) +
  (post_pwp.occ %*% t(propwtr.ub.mean))


marg.wpcp.b.mean <- logistic(as.vector(int.occ) + bwp.b.mean)
marg.wpcp.ub.mean <- logistic(as.vector(int.occ) + bwp.ub.mean)

probs_wpcp.b.mean <- melt(marg.wpcp.b.mean, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.b.mean$burn <- 1

probs_wpcp.ub.mean <- melt(marg.wpcp.ub.mean, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.ub.mean$burn <- 0

probs_wpcp.b.mean$wpcp <- wpcp.b_df$wpcp[probs_wpcp.b.mean$wpid]
probs_wpcp.b.mean$wpcpconv <- (probs_wpcp.b.mean$wpcp * sd(occ.us$medwtpcp)) + mean(occ.us$medwtpcp)

probs_wpcp.ub.mean$wpcp <- wpcp.ub_df$wpcp[probs_wpcp.ub.mean$wpid]
probs_wpcp.ub.mean$wpcpconv <- (probs_wpcp.ub.mean$wpcp * sd(occ.us$medwtpcp)) + mean(occ.us$medwtpcp)

probs_wpcp.mean <- rbind(probs_wpcp.b.mean, probs_wpcp.ub.mean)

prob_wp.occ.mean <- probs_wpcp.mean %>%
  group_by(wpcpconv,burn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )


p_wpcp.occ.mean <- ggplot(prob_wp.occ.mean) + geom_line(aes(y=p, x=wpcpconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=wpcpconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Median winter \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

##Marginal effects of winter precip on occurrence (while holding spring precip at its minimum)
sppcp.occ.b.min = as.matrix(rep(min(occ.all[occ.all$burn== 1,]$medsppcp), length.out=nrow(occ.all[occ.all$burn== 1,])))
sppcp.occ.ub.min = as.matrix(rep(min(occ.all[occ.all$burn== 0,]$medsppcp), length.out=nrow(occ.all[occ.all$burn== 0,])))

#transform back to original scales for proportion
sppcp.occ.b.min.us <- sppcp.occ.b.min * sd(occ.us$medsppcp) + mean(occ.us$medsppcp)
sppcp.occ.ub.min.us <- sppcp.occ.ub.min * sd(occ.us$medsppcp) + mean(occ.us$medsppcp)

#generate proportion and rescale
propwtr.b.min = scale(as.matrix(wtpcp.occ.b.us/(wtpcp.occ.b.us + sppcp.occ.b.min.us)))
propwtr.ub.min = scale(as.matrix(wtpcp.occ.ub.us/(wtpcp.occ.ub.us + sppcp.occ.ub.min.us)))

bwp.b.min <- (post_wp.occ %*% t(wtpcp.occ.b)) + (post_b.occ %*% t(b.occ)) +
  (post_pwp.occ %*% t(propwtr.b.min))
bwp.ub.min <- (post_wp.occ %*% t(wtpcp.occ.ub)) + (post_b.occ %*% t(ub.occ)) +
  (post_pwp.occ %*% t(propwtr.ub.min))

marg.wpcp.b.min <- logistic(as.vector(int.occ) + bwp.b.min)
marg.wpcp.ub.min <- logistic(as.vector(int.occ) + bwp.ub.min)

probs_wpcp.b.min <- melt(marg.wpcp.b.min, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.b.min$burn <- 1

probs_wpcp.ub.min <- melt(marg.wpcp.ub.min, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.ub.min$burn <- 0

probs_wpcp.b.min$wpcp <- wpcp.b_df$wpcp[probs_wpcp.b.min$wpid]
probs_wpcp.b.min$wpcpconv <- (probs_wpcp.b.min$wpcp * sd(occ.us$medwtpcp)) + mean(occ.us$medwtpcp)

probs_wpcp.ub.min$wpcp <- wpcp.ub_df$wpcp[probs_wpcp.ub.min$wpid]
probs_wpcp.ub.min$wpcpconv <- (probs_wpcp.ub.min$wpcp * sd(occ.us$medwtpcp)) + mean(occ.us$medwtpcp)

probs_wpcp.min <- rbind(probs_wpcp.b.min, probs_wpcp.ub.min)

prob_wp.occ.min <- probs_wpcp.min %>%
  group_by(wpcpconv,burn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )

p_wpcp.occ.min <- ggplot(prob_wp.occ.min) + geom_line(aes(y=p, x=wpcpconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=wpcpconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Median winter \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

p_occ.comb <- plot_grid(p_wpcp.occ.h + ylim(c(0,1)),
                        p_wpcp.occ.mean +ylim(c(0,1)) +ylab(""),
                        p_wpcp.occ.min +ylim(c(0,1)) + ylab(""),
                        nrow=1,
                        labels=c("Spring precipitation (max)", "Spring precipitation (mean)", "Spring precipitation (min)"),
                        label_fontfamily = "Times New Roman",
                        label_size = 10,
                        hjust = c(-0.6,-0.6,-0.6),
                        vjust = 1.25)

#Prevalence and precip
#load scaled values
rr.all <- read.csv(paste0(datafolder,"/reprate_model_all_hillshade.csv"), stringsAsFactors = FALSE)

#extract winter precip cov
post_wp.rr <- as.matrix(posterior.rr[,9])
#extract spring precip
post_sp.rr <- as.matrix(posterior.rr[,8])

#extract proportion cov
post_pwp.rr <- as.matrix(posterior.rr[,10])

wtpcp.rr.b = as.matrix(seq(min(rr.all[rr.all$burn== 1,]$wtpr), max(rr.all[rr.all$burn== 1,]$wtpr), length.out=nrow(rr.all[rr.all$burn== 1,])))
wtpcp.rr.ub = as.matrix(seq(min(rr.all[rr.all$burn== 0,]$wtpr), max(rr.all[rr.all$burn== 0,]$wtpr), length.out=nrow(rr.all[rr.all$burn== 0,])))
b.rr = as.matrix(rep(1,nrow(wtpcp.rr.b)))
ub.rr = as.matrix(rep(0,nrow(wtpcp.rr.ub)))
sppcp.rr.b = as.matrix(rep(max(rr.all[rr.all$burn== 1,]$sppr), length.out=nrow(rr.all[rr.all$burn== 1,])))
sppcp.rr.ub = as.matrix(rep(max(rr.all[rr.all$burn== 0,]$sppr), length.out=nrow(rr.all[rr.all$burn== 0,])))

#load unscaled data to get proportion correct
rr.us <- read.csv(paste0(datafolder,"/reprate_model_all_noscale_hillshade.csv"), stringsAsFactors = FALSE)
#generate transformed values
wtpcp.rr.b.us <- wtpcp.rr.b * sd(rr.us$wtpr) + mean(rr.us$wtpr) #convert back to original scale
wtpcp.rr.ub.us <- wtpcp.rr.ub * sd(rr.us$wtpr) + mean(rr.us$wtpr) #convert back to original scale
sppcp.rr.b.us <- sppcp.rr.b * sd(rr.us$sppr) + mean(rr.us$sppr)
sppcp.rr.ub.us <- sppcp.rr.ub * sd(rr.us$sppr) + mean(rr.us$sppr)
#generate proportion and rescale
propwtr.rr.b = scale(as.matrix(wtpcp.rr.b.us/(wtpcp.rr.b.us + sppcp.rr.b.us)))
propwtr.rr.ub = scale(as.matrix(wtpcp.rr.ub.us/(wtpcp.rr.ub.us + sppcp.rr.ub.us)))

bwp.b.rr <- (post_wp.rr %*% t(wtpcp.rr.b)) + (post_b.rr %*% t(b.rr)) +
  (post_pwp.rr %*% t(propwtr.rr.b)) +(post_sp.rr %*% t(sppcp.rr.b))
bwp.ub.rr <- (post_wp.rr %*% t(wtpcp.rr.ub)) + (post_b.rr %*% t(ub.rr)) +
  (post_pwp.rr %*% t(propwtr.rr.ub)) +(post_sp.rr %*% t(sppcp.rr.ub))


marg.wpcp.b.h.rr <- logistic(as.vector(int.rr) + bwp.b.rr)
marg.wpcp.ub.h.rr <- logistic(as.vector(int.rr) + bwp.ub.rr)

probs_wpcp.b.h.rr <- melt(marg.wpcp.b.h.rr, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.b.h.rr$burn <- 1

probs_wpcp.ub.h.rr <- melt(marg.wpcp.ub.h.rr, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.ub.h.rr$burn <- 0

wpcp.b.rr_df <- data.frame(
  wpcp = as.vector(wtpcp.rr.b),
  id = seq(1:length(wtpcp.rr.b))
)

wpcp.ub.rr_df <- data.frame(
  wpcp = as.vector(wtpcp.rr.ub),
  id = seq(1:length(wtpcp.rr.ub))
)

probs_wpcp.b.h.rr$wpcp <- wpcp.b.rr_df$wpcp[probs_wpcp.b.h.rr$wpid]
probs_wpcp.b.h.rr$wpcpconv <- (probs_wpcp.b.h.rr$wpcp * sd(rr.us$wtpr)) + mean(rr.us$wtpr)

probs_wpcp.ub.h.rr$wpcp <- wpcp.ub.rr_df$wpcp[probs_wpcp.ub.h.rr$wpid]
probs_wpcp.ub.h.rr$wpcpconv <- (probs_wpcp.ub.h.rr$wpcp * sd(rr.us$wtpr)) + mean(rr.us$wtpr)

probs_wpcp.h.rr <- rbind(probs_wpcp.b.h.rr, probs_wpcp.ub.h.rr)

prob_wp.rr.high <- probs_wpcp.h.rr %>%
  group_by(wpcpconv,burn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )
##Marginal effects of winter precip on prevalence (while holding spring precip at its maximum)
p_wpcp.rr.h <- ggplot(prob_wp.rr.high) + geom_line(aes(y=p, x=wpcpconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=wpcpconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Annual winter \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

##Marginal effects of winter precip on prevalence (while holding spring precip at its mean)
sppcp.rr.b.mean = as.matrix(rep(mean(rr.all[rr.all$burn== 1,]$sppr), length.out=nrow(rr.all[rr.all$burn== 1,])))
sppcp.rr.ub.mean = as.matrix(rep(mean(rr.all[rr.all$burn== 0,]$sppr), length.out=nrow(rr.all[rr.all$burn== 0,])))

#convert back to unscaled values
sppcp.rr.b.us.mean <- sppcp.rr.b.mean * sd(rr.us$sppr) + mean(rr.us$sppr)
sppcp.rr.ub.us.mean <- sppcp.rr.ub.mean * sd(rr.us$sppr) + mean(rr.us$sppr)
#generate proportion and rescale
propwtr.rr.b.mean = scale(as.matrix(wtpcp.rr.b.us/(wtpcp.rr.b.us + sppcp.rr.b.us.mean)))
propwtr.rr.ub.mean = scale(as.matrix(wtpcp.rr.ub.us/(wtpcp.rr.ub.us + sppcp.rr.ub.us.mean)))

bwp.b.rr.mean <- (post_wp.rr %*% t(wtpcp.rr.b)) + (post_b.rr %*% t(b.rr)) +
  (post_pwp.rr %*% t(propwtr.rr.b.mean)) + (post_sp.rr %*% t(sppcp.rr.b.mean))
bwp.ub.rr.mean <- (post_wp.rr %*% t(wtpcp.rr.ub)) + (post_b.rr %*% t(ub.rr)) +
  (post_pwp.rr %*% t(propwtr.rr.ub.mean)) + (post_sp.rr %*% t(sppcp.rr.ub.mean))


marg.wpcp.b.mean.rr <- logistic(as.vector(int.rr) + bwp.b.rr.mean)
marg.wpcp.ub.mean.rr <- logistic(as.vector(int.rr) + bwp.ub.rr.mean)

probs_wpcp.b.mean.rr <- melt(marg.wpcp.b.mean.rr, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.b.mean.rr$burn <- 1

probs_wpcp.ub.mean.rr <- melt(marg.wpcp.ub.mean.rr, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.ub.mean.rr$burn <- 0

probs_wpcp.b.mean.rr$wpcp <- wpcp.b.rr_df$wpcp[probs_wpcp.b.mean.rr$wpid]
probs_wpcp.b.mean.rr$wpcpconv <- (probs_wpcp.b.mean.rr$wpcp * sd(rr.us$wtpr)) + mean(rr.us$wtpr)

probs_wpcp.ub.mean.rr$wpcp <- wpcp.ub.rr_df$wpcp[probs_wpcp.ub.mean.rr$wpid]
probs_wpcp.ub.mean.rr$wpcpconv <- (probs_wpcp.ub.mean.rr$wpcp * sd(rr.us$wtpr)) + mean(rr.us$wtpr)

probs_wpcp.mean.rr <- rbind(probs_wpcp.b.mean.rr, probs_wpcp.ub.mean.rr)

prob_wp.rr.mean <- probs_wpcp.mean.rr %>%
  group_by(wpcpconv,burn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )

p_wpcp.rr.mean <- ggplot(prob_wp.rr.mean) + geom_line(aes(y=p, x=wpcpconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=wpcpconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Annual winter \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))


##Marginal effects of winter precip on prevalence (while holding spring precip at its minimum)

sppcp.rr.b.min = as.matrix(rep(min(rr.all[rr.all$burn== 1,]$sppr), length.out=nrow(rr.all[rr.all$burn== 1,])))
sppcp.rr.ub.min = as.matrix(rep(min(rr.all[rr.all$burn== 0,]$sppr), length.out=nrow(rr.all[rr.all$burn== 0,])))

#convert back to unscaled values
sppcp.rr.b.us.min <- sppcp.rr.b.min * sd(rr.us$sppr) + mean(rr.us$sppr)
sppcp.rr.ub.us.min <- sppcp.rr.ub.min * sd(rr.us$sppr) + mean(rr.us$sppr)
#generate proportion and rescale
propwtr.rr.b.min = scale(as.matrix(wtpcp.rr.b.us/(wtpcp.rr.b.us + sppcp.rr.b.us.min)))
propwtr.rr.ub.min = scale(as.matrix(wtpcp.rr.ub.us/(wtpcp.rr.ub.us + sppcp.rr.ub.us.min)))

bwp.b.rr.min <- (post_wp.rr %*% t(wtpcp.rr.b)) + (post_b.rr %*% t(b.rr)) +
  (post_pwp.rr %*% t(propwtr.rr.b.min)) +(post_sp.rr %*% t(sppcp.rr.b.min))
bwp.ub.rr.min <- (post_wp.rr %*% t(wtpcp.rr.ub)) + (post_b.rr %*% t(ub.rr)) +
  (post_pwp.rr %*% t(propwtr.rr.ub.min))


marg.wpcp.b.min.rr <- logistic(as.vector(int.rr) + bwp.b.rr.min)
marg.wpcp.ub.min.rr <- logistic(as.vector(int.rr) + bwp.ub.rr.min)

probs_wpcp.b.min.rr <- melt(marg.wpcp.b.min.rr, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.b.min.rr$burn <- 1

probs_wpcp.ub.min.rr <- melt(marg.wpcp.ub.min.rr, value.name = "prob", varnames = c("iter", "wpid"))
probs_wpcp.ub.min.rr$burn <- 0

probs_wpcp.b.min.rr$wpcp <- wpcp.b.rr_df$wpcp[probs_wpcp.b.min.rr$wpid]
probs_wpcp.b.min.rr$wpcpconv <- (probs_wpcp.b.min.rr$wpcp * sd(rr.us$wtpr)) + mean(rr.us$wtpr)

probs_wpcp.ub.min.rr$wpcp <- wpcp.ub.rr_df$wpcp[probs_wpcp.ub.min.rr$wpid]
probs_wpcp.ub.min.rr$wpcpconv <- (probs_wpcp.ub.min.rr$wpcp * sd(rr.us$wtpr)) + mean(rr.us$wtpr)

probs_wpcp.min.rr <- rbind(probs_wpcp.b.min.rr, probs_wpcp.ub.min.rr)

prob_wp.rr.min <- probs_wpcp.min.rr %>%
  group_by(wpcpconv,burn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )

p_wpcp.rr.min <- ggplot(prob_wp.rr.min) + geom_line(aes(y=p, x=wpcpconv, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=wpcpconv, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Annual winter \nprecipitation (mm)") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

p_rr.comb <- plot_grid(p_wpcp.rr.h + ylim(c(0,1)),
                        p_wpcp.rr.mean +ylim(c(0,1)) +ylab(""),
                        p_wpcp.rr.min +ylim(c(0,1)) + ylab(""),
                        nrow=1,
                        labels=c("Spring precipitation (max)", "Spring precipitation (mean)", "Spring precipitation (min)"),
                        label_fontfamily = "Times New Roman",
                        label_size = 10,
                        hjust = c(-0.6,-0.6,-0.6),
                        vjust = 1.25)

precip.plot <- plot_grid(pcp.plot.1  + theme(plot.margin = margin(t=8, unit="pt")),
                         NULL,
                         p_occ.comb,
                         NULL,
                         p_rr.comb,
                         nrow=5,
                         vjust= 1,
                         rel_heights = c(1,0.05,1,0.05,1)
                         )

cowplot::ggsave(paste0(outfolder,"/Figure4.tiff"),plot=precip.plot, width=10.5, height=8, units="in")
