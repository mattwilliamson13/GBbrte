##This code subsets the data for occurrence models for use in models based on all points
library(rethinking) #for coerce_index function
library(tidyverse)

infolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/" #data is already processed so both inputs and outputs in processed folder
outfolder <- "/Users/matthewwilliamson/Google Drive/GB_Final/processed_data/"
# load original data frame
db.or <-read.csv(file =paste0(infolder,"alldatamerged.csv"))
db.or <- db.or[order(db.or$yr),] #order by year

#Subset dataframe for occurrence analysis
db.sub <- db.or[,c(1:5,8:10)]

#subset to determine annual occurrence by point per year
db.sup.p <- db.sub[,3:5]
#reshape into wide format
db.sub.p.r <-reshape(db.sup.p, idvar="pt",timevar = "yr", direction="wide")
#commands to identify first and last occurrence
db.sub.p.r$occind <- apply(db.sub.p.r[,2:14], 1,function(x) substr(colnames(db.sub.p.r)[(min(which(x>0)))+1], start=6, stop=10)) #get the year of first occurrence
db.sub.p.r$absind <- apply(db.sub.p.r[,2:14], 1,function(x) substr(colnames(db.sub.p.r)[(max(which(x == 0)))+1], start=6, stop=10)) #get the year of latest absence
db.sub.p.r$mispres <- db.sub.p.r$occind < db.sub.p.r$absind #if latest absence is greater than first present, then BRTE was detected then lost
db.sub.p.r$lastsamp <- apply(db.sub.p.r[,2:14], 1,function(x) substr(colnames(db.sub.p.r)[(max(which(!is.na(x))))+1], start=6, stop=10)) #get the year of latest absence

#get all sites where BRTE never detected
db.all.abs <- db.sub.p.r[is.na(db.sub.p.r$occind),]
db.all.abs$yr <- db.all.abs$lastsamp #set the year for sites unoccupied throughout the study to the last year measured
db.abs <- db.all.abs[,c(1,19)]


db.pres <- db.sub.p.r[!is.na(db.sub.p.r$occind),]
db.pres$yr <- db.pres$occind #set the year to the earliest occupied year

db.pr <- db.pres[,c(1,19)]
db.occyr <- rbind(db.abs, db.pr)

db.occ.merge <- merge(db.occyr,db.or, by=c("pt","yr")) #merge back to the original dataframe
db.occ.merge$occ <- ifelse(db.occ.merge$brte >= 1, 1,0)

#Calculate proportion of year's precip falling in water
db.occ.merge$propwtr <- db.occ.merge$wtpr/(db.occ.merge$sppr + db.occ.merge$wtpr)


#remove incomplete casess
db.mod <- db.occ.merge[,c(3:4,1,2,17,8:10,15:16,18)]
db.mod.2 <- db.mod[complete.cases(db.mod),] #none of these are correlated in the full dataset

db.occ.merge <- db.mod.2

pt_covs <- read.csv(file=paste0(infolder,"med_pcp_hillshade.csv"))

#Join mean values of precip and aspect
db_join <- db.occ.merge %>% left_join(pt_covs, by=c("range", "area","pt"))
db_join <- db_join[,-c(9:11)] #eliminate "observation level" precip estimates

write.csv(db_join, paste0(outfolder,"occurence_model_all_hillshade_noscale.csv"), row.names=FALSE) #updated for new models

#scale covs
db_join$elev <- scale(db_join$elev)
db_join$medwtpcp <- scale(db_join$medwtpcp)
db_join$medsppcp <- scale(db_join$medsppcp)
db_join$medpropwtr <- scale(db_join$medpropwtr)
db_join$hillshade <- scale(db_join$hillshade)
db.mod.2 <- db_join
#coerce range, canyon, and point to indices
db.mod.2$range_id<- coerce_index(db.mod.2$range)
db.mod.2$area_id<- coerce_index(db.mod.2$area)
db.mod.2$point_id<- coerce_index(db.mod.2$pt)

write.csv(db.mod.2, paste0(outfolder,"occurence_model_all_hillshade.csv"), row.names=FALSE) #updated for new models
