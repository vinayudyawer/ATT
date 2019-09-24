## Croc example

library(tidyverse)
library(sf)
library(raster)
library(mapview)

library(VTrack)
source("~/Documents/GitHub/ATT/R/Functions to be finished/LeastCost UD/2019-08-20_lcDistance.R")

data("PointsCircuitous_crocs")
statinfo <- 
  PointsCircuitous_crocs %>% 
  as_tibble() %>%
  filter(LOCATION > 0) %>% 
  transmute(project_name = "Wenlock",
            installation_name = "Wenlock",
            station_name = LOCATION,
            receiver_name = paste0("VR2W-",LOCATION),
            deploymentdatetime_timestamp = NA,
            recoverydatetime_timestamp = NA,
            status = NA,
            station_longitude = LONGITUDE,
            station_latitude = LATITUDE,
            imos_device = FALSE)

data("crocs")
tagdata <-
  crocs %>%
  as_tibble() %>% 
  transmute(Date.and.Time..UTC. = lubridate::ymd_hms(Date.Time),
            Receiver = Receiver.Name,
            Transmitter = paste(Code.Space, ID, sep="-"),
            Transmitter.Name = Transmitter.Name,
            Transmitter.Serial = Transmitter.S.N,
            Sensor.Value = Sensor.1,
            Sensor.Unit = Units.1,
            Station.Name = Receiver.S.N) %>% 
  left_join(statinfo[c("station_name", "station_longitude","station_latitude")], by=c("Station.Name" = "station_name")) %>% 
  rename(Longitude = station_longitude,
         Latitude = station_latitude)

taginfo <-
  tagdata %>% 
  group_by(Transmitter.Name) %>% 
  slice(1) %>%
  ungroup() %>% 
  transmute(transmitter_id = Transmitter,
            tag_id = c(94, 139, 99),
            tag_project_name = "Wenlock",
            scientific_name = "Crocodylus porosus",
            common_name = "Saltwater Crocodile",
            embargo_date = NA, is_protected = F,
            release_longitude = NA, release_latitude = NA,
            ReleaseDate = NA, sensor_slope = NA, sensor_intercept = NA, sensor_type = NA,
            sensor_unit = NA, tag_model_name = NA, tag_serial_number = Transmitter.Serial,
            tag_expected_life_time_days = 1000, tag_status = "Deployed", sex = NA, measurement = NA,
            dual_sensor_tag = F)

ATTdata <- setupData(Tag.Detections = tagdata %>% filter(Transmitter.Name %in% "Gecko") %>% slice(1:500), 
                     Station.Information = statinfo, 
                     Tag.Metadata = taginfo, 
                     source = "VEMCO")

# abacusPlot(ATTdata, new.window = F)

COAdata <- COA(ATTdata)

## Cost layer

statinfo.sp <-
  statinfo %>% 
  st_as_sf(coords=c("station_longitude", "station_latitude"), crs = 4326) %>% 
  st_transform(crs = 3577)


aus <- st_read(file.choose(), crs = 4326)

map <- 
  aus %>% 
  st_transform(crs = 3577) %>% 
  st_crop(st_bbox(extent(statinfo.sp) + 10000)) %>% 
  st_buffer(dist = 50)

cost.raster <- rasterize(map, raster(extent(map), resolution = 10), 1)
cost.raster[is.na(values(cost.raster))] <- 1000
projection(cost.raster) <- CRS("+init=epsg:3577")

cost <- projectRaster(cost.raster, crs = CRS("+init=epsg:4326"))

plot(cost, col = viridis::viridis(2))
plot(cost.raster, col = viridis::viridis(2))

least.costUD <- lcDistance(ATTdata = ATTdata, ## Station information, tagdata and taginfo data all in one ATTdata object
                           cost = cost,## Cost raster for your study site. If NULL it finds coastline data from OSM server
                           ll_epsg = 4326,    ## EPSG code for the raw data (in lat/long)
                           utm_epsg = 3577,   ## EPSG code for the Projected CRS for your study site (in meters)
                           timestep = 60,     ## Timestep in minutes for COA estimation (see COA() function for details of timestep)
                           h = 100,           ## Smoothing parameter for UD estimate
                           cost.res = 50,     ## Resolution of cost raster used for least cost path estimation
                           UDgrid = 20)       ## Resolution of final UD raster file in meters

plot.lcUD(least.costUD)










