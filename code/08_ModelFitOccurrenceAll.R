#code for actual stan model run
library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
infolder <- here::here("processed_data/") #data is already processed so both inputs and outputs in processed folder
outfolder <- here::here("model_fits/")
stanfolder <- here::here("stan_models/")
db.stan <- read.csv(file=paste0(infolder,"/occurence_model_all_hillshade.csv"), stringsAsFactors = FALSE)
#set up indices for nested effects
can.rg <- db.stan[,13:14]
can.rg.u <- unique(can.rg[,c('area_id','range_id')])
colnames(can.rg.u)[2] <- "CanInRg"
pt.can <- db.stan[,c(14:15,6, 9:12)]
pt.can.u <- unique(pt.can[,c("point_id",'area_id',"elev","hillshade","medwtpcp", "medsppcp", "medpropwtr")])
colnames(pt.can.u)[2] <-c("PtInCan")


## All Data
occ_data <- list(N=nrow(db.stan),R=length(unique(db.stan$range_id)),C=length(unique(db.stan$area_id)),P=length(unique(db.stan$point_id)),
                 PtID=db.stan$point_id, CanID = db.stan$area_id, RgID = db.stan$range_id, CanInRg=can.rg.u$CanInRg, PtInCan=pt.can.u$PtInCan,
                 y=db.stan$occ, W = as.matrix(db.stan[,c(6,9:10,12)]), X = as.matrix(db.stan[,c(7:8)]), J = ncol(db.stan[,c(6,9:10,12)]), K = ncol(db.stan[,c(7:8)]))

occ_fit_all <- stan(paste0(stanfolder,"occurrence_asp.stan"), data=occ_data, seed=082980, control=list(adapt_delta = 0.97))

save(occ_fit_all, file=paste0(outfolder,"occ_fit_all_hillshade.RData"))
