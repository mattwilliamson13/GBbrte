## This code is to generate mean precip values and hillshade values for use as point-level covariates
library(tidyverse)

infolder <- here::here("processed_data/") #data is already processed so both inputs and outputs in processed folder
outfolder <- here::here("processed_data/")
#load the combined dataset
db.or <-read.csv(file =paste0(infolder,"/alldatamerged.csv"))

db.pcp <- db.or %>% group_by(pt) %>% summarise(medwtpcp = median(wtpr, na.rm = TRUE),
                                                medsppcp = median(sppr, na.rm = TRUE))
db.proppcp <- db.or %>% group_by(pt, yr) %>% summarise(propwtr = wtpr/(wtpr+sppr)) %>%
              group_by(pt) %>% summarise(medpropwtr = median(propwtr))

precip <- unique(db.or[,c("range","area", "pt")]) %>%
            left_join(db.pcp, by="pt") %>%
            left_join(db.proppcp, by="pt")
aspect <- read.csv(file=here::here("original_data/Export_vegpts_hillshade_30m.csv"), stringsAsFactors = FALSE) #southwesternnes from Matthias Leu

pt_covs <- aspect[,c(1,2,3,6)] %>% left_join(precip, by=c("range", "area","pt"))
colnames(pt_covs)[4] <- "hillshade"
write.csv(pt_covs, paste0(outfolder,"/med_pcp_hillshade.csv"), row.names=FALSE) #write for use combination with model subset data
