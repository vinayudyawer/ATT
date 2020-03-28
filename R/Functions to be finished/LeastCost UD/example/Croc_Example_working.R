## Croc example

library(tidyverse)
library(sf)
library(raster)
library(mapview)

library(VTrack)
source("2020-03-21_lcDistance.R")

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

Jacko_tagdata <- tagdata %>% filter(Transmitter.Name %in% "Robert D")

ATTdata <- setupData(Tag.Detections = Jacko_tagdata, 
                     Station.Information = statinfo, 
                     Tag.Metadata = taginfo, 
                     source = "VEMCO")

# abacusPlot(ATTdata, new.window = F)
coa_dat <- COA(ATTdata)

coa_sf <-
  coa_dat %>% 
  st_as_sf(coords = c("Longitude.coa", "Latitude.coa"), crs = 4326)

## Cost and Transition layers
wenlock.raster <- raster("data/wenlock raster UTM.tif") %>% ratify()
wenlock.raster[values(wenlock.raster) %in% 0] <- 1000

cost <- projectRaster(wenlock.raster, crs = CRS("+init=epsg:4326"))

trCost <- transition(1/wenlock.raster, mean, directions = 16)
trans <- geoCorrection(trCost, type = "c")

mapview(raster(trans)) + coa_sf

least.costUD <- lcDistance(ATTdata = ATTdata, ## Station information, tagdata and taginfo data all in one ATTdata object
                           trans = trans,     ## Transition layer in UTM (meters)
                           ll_epsg = 4326,    ## EPSG code for the raw data (in lat/long)
                           utm_epsg = 3577,   ## EPSG code for the Projected CRS for your study site (in meters)
                           timestep = 60,     ## Timestep in minutes for COA estimation (see COA() function for details of timestep)
                           h = 100,           ## Smoothing parameter for UD estimate
                           cost.res = 50,     ## Resolution of cost raster used for least cost path estimation
                           UDgrid = 20)       ## Resolution of final UD raster file in meters

plot.lcUD(least.costUD)










