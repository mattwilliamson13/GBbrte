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
occ.ub <- read.csv(paste0(datafolder,"/occurence_model_ub_hillshade_noscale.csv"), stringsAsFactors = FALSE)


occ.graze <- occ.all %>% mutate(dset = "All \npoints")
occ.graze.b <- occ.all %>% filter(., burn == 1) %>% mutate(dset="Burned \npoints")
occ.graze.ub <- occ.all %>% filter(., burn == 0) %>% mutate(dset="Unburned \npoints")
occ.graze <- rbind(occ.graze, occ.graze.b, occ.graze.ub)

#generate a boxplot of actual measured data based on burning and grazing history
graze.comb <-ggplot(data = occ.graze) +
  aes(x = as.factor(graze), fill = as.factor(occ)) +
  geom_bar(stat = "count",colour="black") +
  scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL) +
  scale_x_discrete(name= NULL, labels = c("Not \ngrazed", "Grazed")) +
  facet_wrap( ~ as.factor(dset), strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(y="Number of points") +
  theme(legend.position="none", axis.text.x = element_text(size=10, family="Times New Roman"), text=element_text(size=14,  family="Times New Roman"),
        panel.background = element_blank())

##generate marginal effects for the binary grazing effect
#load occurrence model fits
load(paste0(infolder, "/occ_fit_all_hillshade.RData"))
load(paste0(infolder, "/occ_fit_ub_hillshade.RData"))

#extract the beta posteriors
posterior.occ.all <- as.data.frame(occ_fit_all, pars=c("beta_p", "beta_o"))
posterior.occ.ub <- as.data.frame(occ_fit_ub, pars=c("beta_p", "beta_o"))

#get the intercept
int.occ.all <- as.matrix(occ_fit_all, pars="mu_alpha")
int.occ.ub <- as.matrix(occ_fit_ub, pars="mu_alpha")

#get the burn covariate
post_b.occ.all <- as.matrix(posterior.occ.all[,5])
#no burn covariate for non-burned points

#set up matrices for fitting burn param
burn.occ.all = as.matrix(rep(1,nrow(occ.all[occ.all$burn==1,])))
uburn.occ.all = as.matrix(rep(0,nrow(occ.all[occ.all$burn==0,])))

#create burn effect
bp.occ.all <- post_b.occ.all %*% t(burn.occ.all)
ubp.occ.all <- post_b.occ.all %*% t(uburn.occ.all)

#get grazing cov
post_g.occ.all <- as.matrix(posterior.occ.all[,6])
post_g.occ.ub <- as.matrix(posterior.occ.ub[,5])

##create grazing tables with burn
burn.graze = as.matrix(rep(1, nrow(occ.all[occ.all$burn==1 & occ.all$graze==1,])))
burn.nograze = as.matrix(rep(1, nrow(occ.all[occ.all$burn==1 & occ.all$graze==0,])))

noburn.graze.all = as.matrix(rep(1, nrow(occ.all[occ.all$burn==0 & occ.all$graze==1,])))
noburn.nograze.all = as.matrix(rep(0, nrow(occ.all[occ.all$burn==0 & occ.all$graze==0,])))

noburn.graze.ub = as.matrix(rep(1, nrow(occ.ub[occ.ub$graze==1,])))
noburn.nograze.ub = as.matrix(rep(0, nrow(occ.ub[occ.ub$graze==0,])))

#create graze effect for different burning combinations
bp.burn.graze <- post_g.occ.all %*% t(burn.graze) + post_b.occ.all %*% t(burn.graze)
bp.burn.nograze <- post_b.occ.all %*% t(burn.nograze)
bp.noburn.graze.all <- post_g.occ.all %*% t(noburn.graze.all)
bp.noburn.nograze.all <-  post_g.occ.all %*% t(noburn.nograze.all)

bp.noburn.graze.ub <- post_g.occ.ub %*% t(noburn.graze.ub)
bp.noburn.nograze.ub <- post_g.occ.ub %*% t(noburn.nograze.ub)

#generate marginal posterior predictions
marg.burn.graze <- logistic(as.vector(int.occ.all)  + bp.burn.graze)
marg.burn.nograze <- logistic(as.vector(int.occ.all) + bp.burn.nograze)
marg.noburn.graze.all <- logistic(as.vector(int.occ.all) + bp.noburn.graze.all)
marg.noburn.nograze.all  <- logistic(as.vector(int.occ.all) + bp.noburn.nograze.all)

marg.noburn.graze.ub <- logistic(as.vector(int.occ.ub) + bp.noburn.graze.ub)
marg.noburn.nograze.ub <- logistic(as.vector(int.occ.ub)  + bp.noburn.nograze.ub)

#create dataframes
probs_burn.graze <- melt(marg.burn.graze, value.name = "prob", varnames = c("iter", "id"))
probs_burn.graze$burn <- "Burned"
probs_burn.graze$graze <- "Grazed"
probs_burn.graze$dset <- "All points"

probs_burn.nograze <- melt(marg.burn.nograze, value.name = "prob", varnames = c("iter", "id"))
probs_burn.nograze$burn <- "Burned"
probs_burn.nograze$graze <- "Not grazed"
probs_burn.nograze$dset <- "All points"

probs_noburn.graze.all <- melt(marg.noburn.graze.all, value.name = "prob", varnames = c("iter", "id"))
probs_noburn.graze.all$burn <- "Unburned"
probs_noburn.graze.all$graze <- "Grazed"
probs_noburn.graze.all$dset <- "All points"

probs_noburn.nograze.all <- melt(marg.noburn.nograze.all, value.name = "prob", varnames = c("iter", "id"))
probs_noburn.nograze.all$burn <- "Unburned"
probs_noburn.nograze.all$graze <- "Not grazed"
probs_noburn.nograze.all$dset <- "All points"

probs_noburn.graze.ub <- melt(marg.noburn.graze.ub, value.name = "prob", varnames = c("iter", "id"))
probs_noburn.graze.ub$burn <- "Unburned"
probs_noburn.graze.ub$graze <- "Grazed"
probs_noburn.graze.ub$dset <- "Unburned \npoints"

probs_noburn.nograze.ub <- melt(marg.noburn.nograze.ub, value.name = "prob", varnames = c("iter", "id"))
probs_noburn.nograze.ub$burn <- "Unburned"
probs_noburn.nograze.ub$graze <- "Not grazed"
probs_noburn.nograze.ub$dset <- "Unburned \npoints"

probs.df <- rbind(probs_burn.graze, probs_burn.nograze, probs_noburn.graze.all,
                  probs_noburn.nograze.all, probs_noburn.graze.ub, probs_noburn.nograze.ub)

probs.df$fac <-factor(paste0(probs.df$burn, probs.df$graze), labels=c("Burned/grazed","Burned/not grazed", "Unburned/grazed", "Unburned/not grazed "))



#marginal effects of binary burning and grazing covariates
graze.mar <-ggplot(data = probs.df) +
  aes(x = fac, y=prob, fill=fac) +
  geom_boxplot() +
  #scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL) +
  scale_fill_viridis(discrete=TRUE, option="E", direction= -1, alpha = 0.8)+
  scale_x_discrete(name= NULL) +
  facet_wrap( ~ as.factor(dset), strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(3, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(y="Posterior probability") +
  theme(legend.position="none", axis.text.x = element_text(size=10, family="Times New Roman", angle=90, hjust=1, vjust =0.5), text=element_text(size=14,  family="Times New Roman"),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0, unit = "pt")),
        strip.text.x = element_text(size=14, margin=margin(t=3, b=3)))

graze.plot.1 <- egg::ggarrange(graze.comb, graze.mar,
                               ncol = 2,
                               widths = c(1,1.33))

#load the prevalence data
load(paste0(infolder, "/reprate_fit_all_hillshade.RData"))
load(paste0(infolder, "/reprate_fit_ub_hillshade.RData"))
load(paste0(infolder, "/reprate_fit_b_hillshade2.RData"))

#load the original data frames
db.rr.all <- read.csv(paste0(datafolder,"/reprate_model_all_hillshade.csv"))
db.rr.ub <- read.csv(paste0(datafolder,"/reprate_model_unburned_hillshade.csv"))
db.rr.b <- read.csv(paste0(datafolder,"/reprate_model_burned_hillshade.csv"))

#extract posteriors
posterior.rr.all <- as.data.frame(reprate_fit_all, pars=c("beta_p", "beta_o"))
posterior.rr.ub <- as.data.frame(reprate_fit_ub, pars=c("beta_p", "beta_o"))
posterior.rr.b <- as.data.frame(reprate_fit_b, pars=c("beta_p", "beta_o"))


#get the intercept
int.rr.all <- as.matrix(reprate_fit_all, pars="mu_alpha")
int.rr.ub <- as.matrix(reprate_fit_ub, pars="mu_alpha")
int.rr.b <- as.matrix(reprate_fit_b, pars="mu_alpha")


#extract the burn covariate
post_b.rr <- as.matrix(posterior.rr.all[,5])
#no burn cov for the unburned or burned only points

#set up binary matrices
burn.rr = as.matrix(rep(1,nrow(db.rr.all[db.rr.all$burn==1,])))
uburn.rr = as.matrix(rep(0,nrow(db.rr.all[db.rr.all$burn==0,])))


bp.rr <- post_b.rr %*% t(burn.rr)
ubp.rr <- post_b.rr %*% t(uburn.rr)

#get proportion grazed covariate
post_g.all <- as.matrix(posterior.rr.all[,6])
post_g.ub <- as.matrix(posterior.rr.ub[,5])
post_g.b <- as.matrix(posterior.rr.b[,6])


#set up prediciton data
graze.b.all = as.matrix(seq(min(db.rr.all[db.rr.all$burn == 1,]$pgraze), max(db.rr.all[db.rr.all$burn == 1,]$pgraze), length.out=nrow(db.rr.all[db.rr.all$burn == 1,]))) #mean=.75, sd=.39
graze.ub.all = as.matrix(seq(min(db.rr.all[db.rr.all$burn == 0,]$pgraze), max(db.rr.all[db.rr.all$burn == 0,]$pgraze), length.out=nrow(db.rr.all[db.rr.all$burn == 0,]))) #mean=.75, sd=.39
mean.prop.grz.alldata <- read.csv(here::here("processed_data/reprate_model_all_noscale_hillshade.csv")) %>% 
  summarise(mean = mean(pgraze), stdev = sd(pgraze))


graze.ub <- as.matrix(seq(min(db.rr.ub$pgraze), max(db.rr.ub$pgraze), length.out=nrow(db.rr.ub)))
mean.prop.grz.ub <- read.csv(here::here("processed_data/reprate_model_unburned_noscale_hillshade.csv")) %>% 
  summarise(mean = mean(pgraze), stdev = sd(pgraze))

graze.b <- as.matrix(seq(min(db.rr.b$pgraze), max(db.rr.b$pgraze), length.out=nrow(db.rr.b)))
mean.prop.grz.b <- read.csv(here::here("processed_data/reprate_model_burned_noscale_hillshade.csv")) %>% 
  summarise(mean = mean(pgraze), stdev = sd(pgraze))



bg.b.all <- (post_g.all %*% t(graze.b.all)) + (post_b.rr %*% t(burn.rr))
bg.ub.all <- post_g.all %*% t(graze.ub.all)
bg.ub <- post_g.ub %*% t(graze.ub)
bg.b <- post_g.b %*% t(graze.b)



#inverse link function
graze_b.all <- logistic(as.vector(int.rr.all) + bg.b.all) # burned
graze_ub.all <- logistic(as.vector(int.rr.all) + bg.ub.all)
graze_b <- logistic(as.vector(int.rr.b) + bg.b)
graze_ub <- logistic(as.vector(int.rr.ub) + bg.ub)


#create dataframe
probs_graze_b.all <- melt(graze_b.all, value.name = "prob", varnames = c("iter", "grazeid"))
probs_graze_b.all$burn <- 1
probs_graze_ub.all <- melt(graze_ub.all, value.name = "prob", varnames = c("iter", "grazeid"))
probs_graze_ub.all$burn <- 0

probs_graze_b <- melt(graze_b, value.name = "prob", varnames = c("iter", "grazeid"))
probs_graze_b$burn <- 1
probs_graze_ub <- melt(graze_ub, value.name = "prob", varnames = c("iter", "grazeid"))
probs_graze_ub$burn <- 0




graze.b.all_df <- data.frame(
  graze = as.vector(graze.b.all),
  id = seq(1:length(graze.b.all))
)

graze.ub.all_df <- data.frame(
  graze = as.vector(graze.ub.all),
  id = seq(1:length(graze.ub.all))
)


graze.b_df <- data.frame(
  graze = as.vector(graze.b),
  id = seq(1:length(graze.b))
)

graze.ub_df <- data.frame(
  graze = as.vector(graze.ub),
  id = seq(1:length(graze.ub))
)

# join grazing data to probs_df
probs_graze_b.all$graze <- graze.b.all_df$graze[probs_graze_b.all$grazeid]
probs_graze_ub.all$graze <- graze.ub.all_df$graze[probs_graze_ub.all$grazeid]
probs_df.all <- rbind(probs_graze_b.all, probs_graze_ub.all)

probs_graze_b$graze <- graze.b_df$graze[probs_graze_b$grazeid]
probs_graze_b$propgraze <- as.numeric(mean.prop.grz.b[1,1]) + (probs_graze_b$graze * as.numeric(mean.prop.grz.b[1,2]))

probs_graze_ub$graze <- graze.ub_df$graze[probs_graze_ub$grazeid]
probs_graze_ub$propgraze <- as.numeric(mean.prop.grz.ub[1,1]) + (probs_graze_ub$graze * as.numeric(mean.prop.grz.ub[1,2]))

probs_df.comb <- rbind(probs_graze_b, probs_graze_ub)


# compute summaries and then plot
prob_g.all <- probs_df.all %>%
  group_by(burn,graze) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  ) 
prob_g.all$propgraze <- as.numeric(mean.prop.grz.alldata[1,1]) + (prob_g.all$graze * as.numeric(mean.prop.grz.alldata[1,2]))

prob_g.comb <- probs_df.comb %>%
  group_by(burn,propgraze) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  )


p_graze.all <- ggplot(prob_g.all) + geom_line(aes(y=p, x=propgraze, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=propgraze, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Proportion of years \ngrazed") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))

p_graze.comb <- ggplot(prob_g.comb) + geom_line(aes(y=p, x=propgraze, colour = factor(burn)), lwd=1.5)+
  geom_ribbon(aes(ymin=lo, ymax=hi, x=propgraze, fill = factor(burn)), alpha = 0.2)+
  scale_colour_viridis(discrete = "TRUE", option="E")+
  scale_fill_viridis(discrete = "TRUE", option="E") +
  labs(y="Posterior probability", x = "Proportion of years \ngrazed") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))


graze.plot.2 <- plot_grid(p_graze.all,
                          p_graze.comb,
                          ncol=1,
                          align = "v",
                          axis= "l",
                          labels = c("C", "D"),
                          label_fontfamily = "Times New Roman",
                          label_x = 0.91,
                          label_y= 1,
                          hjust = 0.05,
                          vjust= 1.2)

graze.plot.all <- plot_grid(graze.plot.1,
                            graze.plot.2,
                            nrow=1,
                            rel_widths = c(2.5,1))

graze.plot.all <- graze.plot.all + draw_label("A", x=0.3, y=0.98, size=14, fontfamily = "Times New Roman", fontface="bold") +
  draw_label("B", x=0.68, y=0.98, size=14, fontfamily = "Times New Roman", fontface="bold")

cowplot::ggsave(paste0(outfolder,"/Figure3.tiff"),plot=graze.plot.all, width=10.5, height=6, units="in")
