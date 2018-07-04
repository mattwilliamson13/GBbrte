#This code subsets the data for a binomial model of prevalence. This analysis is aimed at understanding what drives local BRTE prevalence given that it occurs (based on the Occurrence analysis).

library(rethinking) #for coerce_index function
library(tidyverse)
infolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/" #data is already processed so both inputs and outputs in processed folder
outfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/"

# load original data frame
db.or <-read.csv(file =paste0(infolder,"alldatamerged.csv"))
db.or <- db.or[order(db.or$yr),] #order by year

#Subset dataframe for relative rate analysis
db.sub <- db.or[,c(1:6,8:11,13:16)]
db.sub$propwtr <- db.sub$wtpr/(db.sub$sppr + db.sub$wtpr)

#subset to determine annual occurrence by point per year
db.sup.p <- db.sub[,3:5]

db.sub.p.r <-reshape(db.sup.p, idvar="pt",timevar = "yr", direction="wide")
db.sub.p.r$totpres <- rowSums(db.sub.p.r[,2:14], na.rm=TRUE)

#select only points where BRTE was detected at least once (i.e. totpres != 0)
db.pres <- as.character(db.sub.p.r[db.sub.p.r$totpres != 0, 1])

db.sub <- db.sub[db.sub$pt %in% db.pres,]

#select proper columns and remove incomplete cases
db.rr <- db.sub[,c(1:3, 5:8, 10, 12:15)] #drop year, tburn (because it doesn't make sense for unburned sites), and graze (because we are using the propgraze data)
db.rr.cc <- db.rr[complete.cases(db.rr),] #no correlations to be concerned about.

db.sub <- db.rr.cc

pt_covs <- read.csv(file=paste0(infolder,"med_pcp_hillshade.csv"))

#Join mean values of precip and aspect
db_join <- db.sub %>% left_join(pt_covs, by=c("range", "area","pt"))
#remove medsppcp because it is highly correlated with medwtpcp (also done for all occupancy analyses)
db_join <- db_join[,-15]
write.csv(db_join, paste0(outfolder,"reprate_model_all_noscale_hillshade.csv"), row.names=FALSE) #updated for new models

#scale covs
db_join$elev <- scale(db_join$elev)
db_join$wtpr <- scale(db_join$wtpr)
db_join$sppr <- scale(db_join$sppr)
db_join$pgrassf <- scale(db_join$pgrassf)
db_join$pgraze <- scale(db_join$pgraze)
db_join$propwtr <- scale(db_join$propwtr)
db_join$medwtpcp <- scale(db_join$medwtpcp)
db_join$medpropwtr <- scale(db_join$medpropwtr)
db_join$hillshade <- scale(db_join$hillshade)
#create indices usinig rethinking package
db_join$range_id<- coerce_index(db_join$range)
db_join$area_id<- coerce_index(db_join$area)
db_join$point_id<- coerce_index(db_join$pt)
db_join$obs_id <- seq(1,nrow(db_join))

#write out the data
write.csv(db_join, paste0(outfolder, "reprate_model_all_hillshade.csv" ), row.names=FALSE) #updated for new models
