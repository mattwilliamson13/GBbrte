library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
infolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/" #data is already processed so both inputs and outputs in processed folder
outfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/model_fits/"
stanfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/stan_models/"

db.stan <- read.csv(file=paste0(infolder,"reprate_model_burned_hillshade2.csv"), stringsAsFactors = FALSE)

#set up indices for nested effects
can.rg <- db.stan[,18:19]
can.rg.u <- unique(can.rg[,c('area_id','range_id')])
colnames(can.rg.u)[2] <- "CanInRg"
pt.can <- db.stan[,c(19:20,6,14:17)]
pt.can.u <- unique(pt.can[,c("point_id",'area_id',"elev", "hillshade","medwtpcp", "medsppcp","medpropwtr")])
colnames(pt.can.u)[2] <-c("PtInCan")

reprate_data_burned <- list(N=nrow(db.stan),R=length(unique(db.stan$range_id)),C=length(unique(db.stan$area_id)),P=length(unique(db.stan$point_id)),
                PtID=db.stan$point_id, CanID = db.stan$area_id, RgID = db.stan$range_id, CanInRg=can.rg.u$CanInRg, PtInCan=pt.can.u$PtInCan,
                nsamp = db.stan$nsamp, y=db.stan$brte, W = as.matrix(pt.can.u[,c(3:7)]), X = as.matrix(db.stan[,c(7:13)]), J = ncol(pt.can.u[,c(3:7)]), K = ncol(db.stan[,c(7:13)]))

y <- db.stan$brte

reprate_fit_b <- stan(paste0(stanfolder,"reporting_rate_asp.stan"), data=reprate_data_burned, seed=082980, control=list(adapt_delta = 0.999, max_treedepth = 16))

save(reprate_fit_b, file=paste0(outfolder,"reprate_fit_b_hillshade2.RData"))
