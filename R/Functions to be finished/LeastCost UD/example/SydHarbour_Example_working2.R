## Sydney Harbour example
## Working script

library(tidyverse)
library(sf)
library(raster)
library(mapview)
library(gdistance)
library(spatstat)

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

## Run lcDistance function

UDobj <- lcDistance(ATTdata, utm_epsg = 3577)

plot.lcUD(UDobj)












