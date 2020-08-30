## Sydney Harbour example
## Working script

library(tidyverse)
library(sf)
library(raster)
library(mapview)
library(gdistance)

## Installing VTrack from github if you dont have the recent version
# devtools::install_github("rossdwyer/VTrack")
library(VTrack)
source("2020-08-29_lcDistance.R")

# Input example datasets
map <- readRDS("data/Example_map.RDS")
statinfo <- readRDS("data/Example_statinfo.RDS")
tagdata <- readRDS("data/Example_tagdata.RDS")
taginfo <- readRDS("data/Example_taginfo.RDS")

ATTdata <- setupData(Tag.Detections = tagdata, 
                     Tag.Metadata = taginfo,
                     Station.Information =  statinfo, 
                     source = "IMOS")


### Convert ATT to residence
residence <- ATT_to_residence(ATTdata = ATTdata, 
                              .iResidenceThreshold = 1, 
                              .iTimeThreshold = 3600, 
                              .iCores = 4)


### Create line segments of non-residence events
statinfo <- ATTdata$Station.Information %>% st_as_sf(coords = c("Station.Longitude", "Station.Latitude"), crs = 4326, remove = F)

costras <- extract_costras(statinfo = statinfo, 
                           utm_epsg = 3577, 
                           cost.res = 50)

trans <- cost_to_trans(cost = costras, 
                       utm_epsg = 3577, 
                       cost.res = 50, 
                       directions = 16) 

traj <- nonresidence_to_traj(residence = residence, 
                             trans = trans, 
                             utm_epsg = 3577, 
                             ll_epsg = 4326)






