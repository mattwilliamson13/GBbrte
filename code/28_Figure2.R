library(extrafont)
library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(rethinking)
library(egg)
library(viridis)
infolder <- here::here("model_fits/")
outfolder <- here::here("summary_plots/")
datafolder <- here::here("processed_data/")


####Generate marginal effects plots for burned (binary) predictor
#load the model fits
load(paste0(infolder, "/reprate_fit_all_hillshade.RData"))
load(paste0(infolder, "/occ_fit_all_hillshade.RData"))

#load the original data frames
db.rr <- read.csv(paste0(datafolder,"/reprate_model_all_hillshade.csv"))
db.occ <- read.csv(paste0(datafolder,"/occurence_model_all_hillshade.csv"))

posterior.rr <- as.data.frame(reprate_fit_all, pars=c("beta_p", "beta_o"))
posterior.occ <- as.data.frame(occ_fit_all, pars=c("beta_p", "beta_o"))

#get the intercept
int.rr <- as.matrix(reprate_fit_all, pars="mu_alpha")
int.occ <- as.matrix(occ_fit_all, pars="mu_alpha")

#extract the binary covariate
post_b.rr <- as.matrix(posterior.rr[,5])
post_b.occ <- as.matrix(posterior.occ[,5])

#set up binary matrices
burn.rr = as.matrix(rep(1,nrow(db.rr[db.rr$burn==1,])))
uburn.rr = as.matrix(rep(0,nrow(db.rr[db.rr$burn==0,])))
burn.occ = as.matrix(rep(1,nrow(db.occ[db.occ$burn==1,])))
uburn.occ = as.matrix(rep(0,nrow(db.occ[db.occ$burn==0,])))



bp.rr <- post_b.rr %*% t(burn.rr)
ubp.rr <- post_b.rr %*% t(uburn.rr)
bp.occ <- post_b.occ %*% t(burn.occ)
ubp.occ <- post_b.occ %*% t(uburn.occ)

#generate marginal posterior predictions
burn.eff.rr <- logistic(as.vector(int.rr) +  bp.rr)
uburn.eff.rr <- logistic(as.vector(int.rr) +  ubp.rr)
burn.eff.occ <- logistic(as.vector(int.occ) +  bp.occ)
uburn.eff.occ <- logistic(as.vector(int.occ) +  ubp.occ)

#create dataframe
probs_b.rr <- melt(burn.eff.rr, value.name = "prob", varnames = c("iter", "burnid"))
probs_b.rr$burn <- 1
probs_ub.rr <- melt(uburn.eff.rr, value.name = "prob", varnames = c("iter", "burnid"))
probs_ub.rr$burn <- 0

probs_b.occ <- melt(burn.eff.occ, value.name = "prob", varnames = c("iter", "burnid"))
probs_b.occ$burn <- 1
probs_ub.occ <- melt(uburn.eff.occ, value.name = "prob", varnames = c("iter", "burnid"))
probs_ub.occ$burn <- 0

probs_burn.rr <- rbind(probs_b.rr, probs_ub.rr)
probs_burn.rr$dset <- "Prevalence"
probs_burn.occ <- rbind(probs_b.occ, probs_ub.occ)
probs_burn.occ$dset <- "Occurrence"
probs_burn <- rbind(probs_burn.rr, probs_burn.occ)

burn.mar <-ggplot(data = probs_burn) +
  aes(x = as.factor(burn), y=prob, fill=as.factor(burn)) +
  geom_boxplot() +
  #scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL) +
  scale_fill_viridis(discrete=TRUE, option="E", alpha=0.8)+
  scale_x_discrete(name= NULL, labels = c("Unburned", "Burned")) +
  facet_wrap( ~ as.factor(dset), strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  labs(y="Posterior probability") +
  theme(legend.position="none", axis.text.x = element_text(size=10, family="Times New Roman"), text=element_text(size=14,  family="Times New Roman"),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0, unit = "pt")))


occ.all <- read.csv(paste0(datafolder,"/occurence_model_all_hillshade_noscale.csv"), stringsAsFactors = FALSE)
occ.all$status <- factor(occ.all$burn, labels=c("Not burned points", "Burned points"))

#generate boxplot of observed data based on occurrence
burn.comb <- ggplot() +
  geom_bar(data=occ.all, aes(x=as.factor(status), fill=as.factor(occ)),colour="black", stat="count", position = "dodge") +
  scale_fill_manual(values = c("white", "black"),name=NULL, labels=NULL)+
  scale_x_discrete(name = NULL, labels = c("Unburned \npoints","Burned \npoints")) +
  labs(y="Number of points") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"),
        panel.background = element_blank()) 
  
#combine with marginal effects plot
burn.plot.1 <- egg::ggarrange(burn.comb, burn.mar,
                              ncol = 2
                              )
###Generate marginal effects for time since burn param
load(paste0(infolder, "/reprate_fit_b_hillshade2.RData")) #load model fit data
rr.burned <- read.csv(paste0(datafolder,"/reprate_model_burned_hillshade.csv"))
rr.burned.noscale <- read.csv(paste0(datafolder,"/reprate_model_burned_noscale_hillshade.csv"))
mean.tb.orig <- mean(rr.burned.noscale$tburn)
sd.tb.orig <- sd(rr.burned.noscale$tburn)
#Extract posterior samples for each slope parameter
posterior.tb <- as.data.frame(reprate_fit_b, pars=c("beta_p", "beta_o"))

#plotting "marginal" effects
#get the intercept
int.tb <- as.matrix(reprate_fit_b, pars="mu_alpha")

post_tb <- as.matrix(posterior.tb[,7])
post_tb2 <- as.matrix(posterior.tb[,8])

#set up prediciton data
tburn = as.matrix(seq(min(rr.burned$tburn), max(rr.burned$tburn), length.out=nrow(rr.burned))) 

btb <- post_tb %*% t(tburn)
btb2 <- post_tb2 %*% t(tburn)

#inverse link function
p_tb <- logistic(as.vector(int.tb) + btb + btb2)

#create dataframe
probs_df<- melt(p_tb, value.name = "prob", varnames = c("iter", "tbid"))

#create dataframe of elev for binding
tb_df <- data.frame(
  tburn = as.vector(tburn),
  id = seq(1:length(tburn))
)

# join graze data to probs_df
probs_df$tburn <- tb_df$tburn[probs_df$tbid]

# compute summaries and then plot
prob_tb <- probs_df %>%
  group_by(tburn) %>%
  summarise(
    p = mean(prob),
    lo = quantile(prob, probs=0.1),
    hi = quantile(prob, probs=0.9)
  ) %>% 
  mutate(burnyr = mean.tb.orig + (.$tburn*sd.tb.orig))

p_tburn.b <- ggplot(prob_tb) + geom_line(aes(y=p, x=burnyr), color = cividis(1, begin = 1, alpha=0.8) , lwd=1.5) +
  geom_ribbon(aes(ymin=lo, ymax=hi, x=burnyr), fill = cividis(1, begin=1), alpha = 0.3) +
  scale_x_continuous(limits = c(0,15))+
  labs(y="Posterior probability", x = "Years since fire") +
  theme(legend.position="none", text=element_text(size=14,  family="Times New Roman"), panel.background = element_rect(fill = "white", colour = "grey50"))


burn.plot.all <- plot_grid(burn.plot.1,
                           p_tburn.b + ylim(0,1),
                           align= "h",
                           axis="b",
                           rel_widths = c(2, 1))

burn.plot.all <- burn.plot.all + draw_label("A", x=0.32, y=0.98, size=14, fontfamily = "Times New Roman", fontface="bold") +
  draw_label("B", x=0.64, y=0.98, size=14, fontfamily = "Times New Roman", fontface="bold") +
  draw_label("C", x= 0.97, y=0.98, size=14, fontfamily = "Times New Roman", fontface="bold")

cowplot::ggsave(paste0(outfolder,"/Figure2.tiff"),plot=burn.plot.all, width=10.5, height=6, units="in")
