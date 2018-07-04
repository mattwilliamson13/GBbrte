library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
infolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/" #data is already processed so both inputs and outputs in processed folder
outfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/model_fits/"
stanfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/stan_models/"
db.stan <- read.csv(file=paste0(infolder,"occurence_model_unburned_hillshade.csv"), stringsAsFactors = FALSE)

#set up indices for nested effects
can.rg <- db.stan[,12:13]
can.rg.u <- unique(can.rg[,c('area_id','range_id')])
colnames(can.rg.u)[2] <- "CanInRg"

pt.can <- db.stan[,c(13:14, 6, 8:11)]
pt.can.u <- unique(pt.can[,c("point_id",'area_id',"elev","hillshade","medwtpcp", "medsppcp", "medpropwtr")])
colnames(pt.can.u)[2] <-c("PtInCan")


occ_data_ub <- list(N=nrow(db.stan),R=length(unique(db.stan$range_id)),C=length(unique(db.stan$area_id)),P=length(unique(db.stan$point_id)),
            PtID=db.stan$point_id, CanID = db.stan$area_id, RgID = db.stan$range_id, CanInRg=can.rg.u$CanInRg, PtInCan=pt.can.u$PtInCan,
            y=db.stan$occ, W = as.matrix(db.stan[,c(6,8:9,11)]), X = db.stan[,7], J = ncol(db.stan[,c(6,8:9,11)]), K = 1) #get rid of spring precip again


occ_fit_ub <- stan(paste0(stanfolder,"occurrence_asp_ub.stan"), data=occ_data_ub, seed=082980, control=list(adapt_delta = 0.95))

save(occ_fit_ub, file=paste0(outfolder,"occ_fit_ub_hillshade.RData"))
