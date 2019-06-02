#This code subsets the data for a binomial model of reporting rate. This analysis is aimed at understanding what drives local BRTE prevalence given that it occurs (based on the Occurrence analysis) on burned sites.
library(tidyverse)
library(rethinking) #for coerce_index function
infolder <- here::here("processed_data/") #data is already processed so both inputs and outputs in processed folder
outfolder <- here::here("processed_data/")

# load original data frame
db.or <-read.csv(file =paste0(infolder,"/alldatamerged.csv"))
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

#select proper columns and remove incomplete casess
db.rr.b <- db.sub[db.sub$burn == 1,]

#select proper columns and remove incomplete cases
db.rr.b <- db.rr.b[,c(1:3, 5:7, 10, 11:15)] #drop year and graze (because we are using the propgraze data)

db.rr.cc <- db.rr.b[complete.cases(db.rr.b),] #no major correlations, retaining all covs
db.rr.cc$tburn2 <- db.rr.cc$tburn * db.rr.cc$tburn
db.sub <- db.rr.cc

pt_covs <- read.csv(file=paste0(infolder,"/med_pcp_hillshade.csv"))

#Join mean values of precip and aspect
db_join <- db.sub %>% left_join(pt_covs, by=c("range", "area","pt")) #no correlations

write.csv(db_join, paste0(outfolder, "/reprate_model_burned_noscale_hillshade.csv" ), row.names=FALSE)
#scale covs
db_join$elev <- scale(db_join$elev)
db_join$wtpr <- scale(db_join$wtpr)
db_join$pgrassf <- scale(db_join$pgrassf)
db_join$pgraze <- scale(db_join$pgraze)
db_join$propwtr <- scale(db_join$propwtr)
db_join$tburn <- scale(db_join$tburn)
db_join$tburn2 <- scale(db_join$tburn2)
db_join$sppr <- scale(db_join$sppr)
db_join$medwtpcp <- scale(db_join$medwtpcp)
db_join$medsppcp <- scale(db_join$medsppcp)
db_join$medpropwtr <- scale(db_join$medpropwtr)
db_join$hillshade <- scale(db_join$hillshade)
#re-order so quadratic term is adjacent to tburn term for ease of interp
db_join <- db_join[,c(1:8,13,9:12,14:17)]

#scale covs usinig rethinking package
db_join$range_id<- coerce_index(db_join$range)
db_join$area_id<- coerce_index(db_join$area)
db_join$point_id<- coerce_index(db_join$pt)
db_join$obs_id <- seq(1,nrow(db_join))

#write out the data
write.csv(db_join, paste0(outfolder, "/reprate_model_burned_hillshade.csv" ), row.names=FALSE) #updated for new models
