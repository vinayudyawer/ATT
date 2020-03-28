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
source("2020-03-28_lcDistance.R")

# Input example datasets
map <- readRDS("data/Example_map.RDS")
statinfo <- readRDS("data/Example_statinfo.RDS")
tagdata <- readRDS("data/Example_tagdata.RDS")
taginfo <- readRDS("data/Example_taginfo.RDS")

ATTdata <- setupData(Tag.Detections = tagdata, 
                     Tag.Metadata = taginfo,
                     Station.Information =  statinfo, 
                     source = "IMOS")

# Quick visualisation of the raw data
statinfo.sp <- 
  statinfo %>% 
  st_as_sf(coords=c("station_longitude", "station_latitude"), crs=4326)

coa.data <-
  COA(ATTdata) %>% 
  st_as_sf(coords=c("Longitude.coa", "Latitude.coa"), crs=4326)

mapview(statinfo.sp, layer = "Receivers", legend = F, homebutton = F, alpha.regions = 0) +
  mapview(map, alpha = 0, col.regions = "grey", legend =F, layer = "Coastline polygon", homebutton = F) +
  mapview(coa.data, layer = "COA positions", legend =F, homebutton = F, alpha = 0, col.regions = 2, size = 2)

## If your using your own cost layer, this is how to make one from a polygon shapefile of land/coastline
cost.raster <- rasterize(map, raster(extent(map), resolution = 0.001), 1000)
cost.raster[is.na(values(cost.raster))] <- 1
projection(cost.raster) <- CRS("+init=epsg:4326")

cost.raster %>% 
  rasterToPoints(.) %>%
  as_tibble() %>% 
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = factor(layer))) +
  coord_fixed(expand = 0) +
  scale_fill_discrete(name = "Cost Value", labels = c("1 (Water)", "1000 (Land)")) +
  theme_bw()

## Create a transition layer to speed up UD estimation if using lots of individuals in the same study site
cost.in_utm <- projectRaster(cost.raster, crs = CRS("+init=epsg:3577"), method = "ngb")
cost.ras <- resample(cost.in_utm, raster(extent(cost.in_utm), res = 50), method = "ngb")
projection(cost.ras) <- CRS("+init=epsg:3577")

trCost <- transition(1/cost.ras, mean, directions = 16)
trans.utm <- geoCorrection(trCost, type = "c")

trans.utm %>% 
  raster() %>% 
  rasterToPoints(.) %>%
  as_tibble() %>% 
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = layer)) +
  coord_fixed(expand = 0) +
  theme_bw()


## least-cost UD estimation
least.costUD <- lcDistance(ATTdata = ATTdata,  ## Station information, tagdata and taginfo data all in one ATTdata object
                           trans = trans.utm,  ## Transition layer, function will use this over cost raster if provided to save time
                           cost = cost.raster, ## Cost raster for your study site. If NULL it finds coastline data from OSM server
                           ll_epsg = 4326,     ## EPSG code for the raw data (in lat/long)
                           utm_epsg = 3577,    ## EPSG code for the Projected CRS for your study site (in meters)
                           timestep = 60,      ## Timestep in minutes for COA estimation (see COA() function for details of timestep)
                           h = 100,            ## Smoothing parameter for UD estimate
                           cost.res = 50,      ## Resolution of cost raster used for least cost path estimation
                           UDgrid = 20)        ## Resolution of final UD raster file in meters


summary(least.costUD)

## Raw tagdata with least cost distance calculated between each detection (in meters)
least.costUD$tagdata

## Center of Activity data calculated in the function
mapview(least.costUD$coa.data)

## Kernel areas in square meters for 50% and 95% UDs
least.costUD$kernel.areas

## Spatial object for the least-cost trajectories calculated for each consecutive detection
least.costUD$lc.traj

## UD raster output which can be used to plot contours or maps with
mapview(least.costUD$UD.raster)

## cost raster used to calculate least-cost paths
mapview(least.costUD$cost.raster)

## Quickly plotting the result
plot(least.costUD$UD.raster, col = viridis::viridis(10), zlim=c(0,100))
contour(least.costUD$UD.raster, levels = c(50, 95), add=T)
plot(map, add=T, border = grey(0.5))
plot(st_geometry(statinfo.sp), pch = 20, col="red", cex=0.5, add=T)

## Visualise least-cost UD
plot.lcUD(least.costUD)


