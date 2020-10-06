#### Least cost path calculation using passive telemetry datasets
## Written by VU 2018-02-12
## Last updated: 2020-08-29 (VU)
##
## Function broken into multiple smaller functions to calculate least cost paths between detections from an ATTdata object (from VTrack)
## around landmasses. Uses 'osmdata' to construct cost raster if none is provided
## Calculates residence events and use residence logs to speed up least cost estimates
## ************************************************************************************************************************************* ##
lcDistance <- function(ATTdata, cost = NULL, trans = NULL, utm_epsg, ll_epsg = 4326, timestep = 60, 
                       h = 200, UDextent = 1, UDgrid = 50, cost.res = 50, directions = 16, 
                       residence_threshold = 1, time_threshold = 3600, cores = 2, stepwise = F, ...){
  ####################################################################################
  ## ATTdata        ATTdata object that consists of Station.Information, Tag.Detections and
  ##                Tag.Metadata (output from setupData function in VTrack)
  ## cost           cost raster, with land pixel values of 1000 and sea as 1, if none provided
  ##                this is extracted and calculated using polygon downloaded from 'osmdata' package
  ## trans          Transition layer used to calculate least cost paths. Providing this will make the estimation faster for
  ##                larger datasets. Bypasses transition matrix calculation every time.
  ## utm_epsg       EPSG code for CRS object with projected coordinate reference system (units m)
  ## ll_epsg        EPSG code for CRS object with unprojected coordinate reference system (units deg) [default 4326]
  ## timestep       timestep for Center of activity positions (see COA() function in VTrack)
  ## h              smoothing parameter for kernel density estimator in meters (default 100m)
  ## UDextent       a value controlling the extent of the grid used for the estimation of kernel density (default = 4)
  ## UDgrid         a number giving the size of the grid on which the kernel density should be estimated (in meters; default 100)
  ## cost.res       resolution of cost layer (in degrees) calculated if not provided [default 100]
  ## directions     number of directional axes to consider when calculating
  ##                least cost path (options: 4, 8, 16 [default])
  ####################################################################################
  
  ## load required libraries and set up CRSs
  sapply(c("lubridate","sf","dplyr","raster","gdistance","spatstat","maptools","rgeos","VTrack", 
           "data.table", "lwgeom"), require, character.only=TRUE, warn.conflicts = F, quietly = T)
  
  ll <- CRS(paste0("+init=epsg:", ll_epsg))
  utm <- CRS(paste0("+init=epsg:", utm_epsg))
  
  ### Convert ATT to residence
  message("- Constructing residence log and determining nonresidence events")
  residence <- ATT_to_residence(ATTdata = ATTdata, 
                                .iResidenceThreshold = residence_threshold, 
                                .iTimeThreshold =  time_threshold,
                                .iCores = cores)
  
  ### Create spatial layers of statinfo and detection rates
  statinfo <- 
    ATTdata$Station.Information %>% 
    st_as_sf(coords = c("Station.Longitude", "Station.Latitude"), crs = 4326, remove = F)
  
  ### Configure Cost or Transition layer
  if(is.null(trans)){
    if(is.null(cost)){
      message("- No cost or transition layer provided. Downloading coastline data from Open Street Map server")
      suppressWarnings(
        cost <- extract_costras(statinfo = statinfo,
                                utm_epsg = utm_epsg, 
                                cost.res = cost.res) 
      )
      suppressWarnings(
        trans <- cost_to_trans(cost = cost, 
                               utm_epsg = utm_epsg, 
                               cost.res = cost.res, 
                               directions = directions)  
      )
    } else {
      message("- No transition layer provided. Accessing cost layer provided")
      if(!class(cost) %in% "RasterLayer"){
        stop("Make sure input cost layer is a Raster object")
      }
      cost.in_utm <- suppressWarnings(raster::projectRaster(cost, crs = utm, method = "ngb"))
      cost <- resample(cost.in_utm, raster(extent(cost.in_utm), res = cost.res), method = "ngb")
      projection(cost) <- utm
      
      suppressWarnings(
        trans <- cost_to_trans(cost = cost, 
                               utm_epsg = utm_epsg, 
                               cost.res = cost.res, 
                               directions = directions)  
      )
    }
  } else {
    message("- Accessing provided transition layer for least cost path estimation")
    if(!class(trans) %in% "TransitionLayer"){
      stop("Make sure input transition layer is a Transition Layer")
    }
    trans <- trans
  }
  
  ## Construct shortest path trajectories between sequence of detection steps
  traj <- nonresidence_to_traj(residence = residence, 
                               trans = trans, 
                               utm_epsg = utm_epsg, 
                               ll_epsg = ll_epsg)
  
  ## Convert trajectories to points
  # pts <- traj_to_points(traj)
  
  ## UD estimation
  suppressWarnings(
    UDras <- ud_est(ATTdata = ATTdata, 
                    traj = traj, 
                    residence = residence, 
                    cost.ras = cost, 
                    ll_epsg = ll_epsg, 
                    utm_epsg = utm_epsg, 
                    timestep = timestep, 
                    UDextent = UDextent, 
                    UDgrid = UDgrid, 
                    h = h, 
                    stepwise = stepwise) 
  )
  
  ## Prep output data
  tagdata <-
    residence$residenceslog %>% 
    left_join(ATTdata$Tag.Metadata[c("Transmitter", "Tag.ID", "Sci.Name", "Common.Name", "Sex", "Bio")], 
              by = c("TRANSMITTERID" = "Transmitter")) %>% 
    transmute(date_time = DATETIME,
              residence_event = RESIDENCEEVENT,
              transmitter_id = TRANSMITTERID,
              station_name = STATIONNAME,
              residence_duration_sec = ELAPSED,
              longitude = Station.Longitude,
              latitude = Station.Latitude,
              tag_id = Tag.ID,
              sci_name = Sci.Name,
              common_name = Common.Name,
              sex = Sex, bio = Bio)
    
  
  out <- list(tagdata = tagdata,
              kernel.areas = data.frame(UD50_m2 = UDras$UD50_area, UD95_m2 = UDras$UD95_area),
              lc.traj = traj %>% st_transform(crs = ll_epsg),
              UD.raster = UDras$raster_output)
  
  return(out)
  
}


## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
## Function to construct residence log from ATT object
ATT_to_residence <- function(ATTdata, .iResidenceThreshold = 2, 
                             .iTimeThreshold = 3600, .iCores = 4){
  
  tags_to_include <- unique(ATTdata$Tag.Metadata$Transmitter)
  
  teldat <-
    ATTdata$Tag.Detections %>%
    transmute(DATETIME = Date.Time,
              TRANSMITTERID = Transmitter,
              SENSOR1 = Sensor.Value,
              UNITS1 = Sensor.Unit,
              RECEIVERID = Receiver,
              STATIONNAME = Station.Name) %>%
    filter(TRANSMITTERID %in% tags_to_include) %>% 
    data.frame()
  
  residence <- RunResidenceExtraction(sInputFile = teldat,
                                      sLocation = "STATIONNAME",
                                      iResidenceThreshold = .iResidenceThreshold,
                                      iTimeThreshold = .iTimeThreshold, 
                                      iCores = .iCores)
  
  res_out <-
    residence %>%
    map(as_tibble)
  
  res_out$residences <-
    res_out$residences %>%
    arrange(STARTTIME) %>% 
    left_join(ATTdata$Station.Information[c("Station.Name", "Station.Longitude", "Station.Latitude")], 
              by = c("STATIONNAME" = "Station.Name"))
  
  res_out$residenceslog <-
    res_out$residenceslog %>%
    arrange(DATETIME) %>% 
    left_join(ATTdata$Station.Information[c("Station.Name", "Station.Longitude", "Station.Latitude")], 
              by = c("STATIONNAME" = "Station.Name"))
  
  res_out$nonresidences <-
    res_out$nonresidences %>%
    arrange(STARTTIME) %>% 
    left_join(ATTdata$Station.Information[c("Station.Name", "Station.Longitude", "Station.Latitude")], 
              by = c("STATIONNAME1" = "Station.Name")) %>% 
    rename(Start.Longitude = Station.Longitude,
           Start.Latitude = Station.Latitude)  %>% 
    left_join(ATTdata$Station.Information[c("Station.Name", "Station.Longitude", "Station.Latitude")], 
              by = c("STATIONNAME2" = "Station.Name")) %>% 
    rename(End.Longitude = Station.Longitude,
           End.Latitude = Station.Latitude)
  
  return(res_out)
}
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
## Extract landmass to calculate transition and cost raster if not provided
extract_costras <- function(statinfo, utm_epsg, cost.res = 100) {
  
  require("osmdata")
  tryCatch({
    
    utm <- CRS(paste0("+init=epsg:", utm_epsg))
    
    bb <- st_bbox(extent(statinfo) + 0.1)
    
    polydat <-
      opq(bbox = bb) %>%
      add_osm_feature(key = 'natural', value = c('water', 'coastline')) %>%
      osmdata_sf(quiet = T)
    
    islands <-
      polydat$osm_polygons %>%
      filter(!natural %in% 'water') %>% 
      rowid_to_column() %>%
      mutate(area = st_area(geometry)) %>%
      dplyr::select(rowid, area) %>% 
      st_combine() %>% st_sf()
    
    water <- 
      polydat$osm_lines %>% 
      st_union() %>% st_line_merge() %>% 
      st_polygonize() %>% st_collection_extract() %>% st_union() %>% st_sf()
    
    coastline <- polydat$osm_lines %>% filter(natural %in% 'coastline') %>% st_union %>% st_line_merge
    pol <- st_as_sfc(bb) %>% st_sf(crs = 4326)
    seapoly <- lwgeom::st_split(st_geometry(pol), st_geometry(coastline)) %>% st_collection_extract() %>% st_sf()
    
    sea <-
      seapoly %>% 
      rowid_to_column() %>%
      mutate(area = st_area(geometry)) %>%
      filter(!area %in% max(area)) %>% 
      bind_rows(water) %>% 
      st_union()
    
    studypol <- 
      sea %>% 
      st_difference(islands) %>%
      st_crop(pol)
    
    poly <-
      studypol %>%
      st_transform(utm_epsg) %>%
      st_crop(extent(statinfo %>% st_transform(utm_epsg)) + 5000) %>% 
      st_sf()
    
    cost.ras <- raster::rasterize(poly, 
                                  raster(extent(statinfo %>% st_transform(utm_epsg)) + 5000, 
                                         res = cost.res), 1)
    cost.ras[is.na(values(cost.ras))] <- 1000
    projection(cost.ras) <- utm
      
    }, error=function(e){message("Error: no land in sight!\nConsider adding your own cost layer")})
  
  return(cost.ras)
}
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
## Convert cost layer to transition layer if required
cost_to_trans <- function(cost, utm_epsg, cost.res = 100, directions = 16){
  
  utm <- CRS(paste0("+init=epsg:", utm_epsg))
  
  if(isLonLat(cost)){
    message("- Projecting input cost layer to UTM")
    cost.in_utm <- raster::projectRaster(cost, crs = utm, method = "ngb")
  } else {
    cost.in_utm <- cost
  }
  
  cost.ras <- resample(cost.in_utm, raster(extent(cost.in_utm), res = cost.res), method = "ngb")
  projection(cost.ras) <- utm
  
  ## Produce transition matrices, and correct for distortion
  tryCatch({
    message("- Calculating transition matrices for least cost path estimation")
    trCost <- transition(1/cost.ras, mean, directions = directions)
    trCost <- geoCorrection(trCost, type = "c")
  }, error=function(e){message("Error in calculating Transition layer")})
  
  return(trCost)
}
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
## Function to generate least cost trajectories from non-residence logs
nonresidence_to_traj <- function(residence, trans, utm_epsg, ll_epsg = 4326){
  
  ll <- CRS(paste0("+init=epsg:", ll_epsg))
  utm <- CRS(paste0("+init=epsg:", utm_epsg))
  
  dat_sf <-
    residence$residenceslog %>% 
    st_as_sf(coords = c("Station.Longitude", "Station.Latitude"), 
             crs = ll_epsg) %>% 
    st_transform(crs = utm_epsg)
  
  ## Setup input data for shortest path trajectory estimation
  ## Select unique combinations of movements to reduce time
  traj_df <- 
    residence$nonresidences %>% 
    dplyr::distinct(STATIONNAME1, STATIONNAME2, .keep_all = T) %>% 
    dplyr::select(-c(STARTTIME, ENDTIME, NONRESIDENCEEVENT, 
                     TRANSMITTERID, DURATION, DISTANCE, ROM))
  
  
  ## Construct shortest path trajectories between each unique combination
  tryCatch({
    message("- Constructing least cost path trajectories between consecutive detections")
    for(i in 1:nrow(traj_df)){
      if(i %in% 1){pb <- txtProgressBar(min = 1, max = nrow(traj_df), style = 3)}
      
      origin <- 
        traj_df %>% slice(i) %>% 
        st_as_sf(coords = c("Start.Longitude", "Start.Latitude"), 
                 crs = ll_epsg, remove = F) %>% 
        st_transform(crs = utm_epsg) %>% 
        as_Spatial()
      
      goal <- 
        traj_df %>% slice(i) %>% 
        st_as_sf(coords = c("End.Longitude", "End.Latitude"), 
                 crs = ll_epsg, remove = F) %>% 
        st_transform(crs = utm_epsg) %>% 
        as_Spatial()
      
      if(i %in% 1){
        paths <-
          shortestPath(trans, origin, goal, output = "SpatialLines") %>% 
          st_as_sf() %>% 
          mutate(STATIONNAME1 = origin$STATIONNAME1,
                 STATIONNAME2 = origin$STATIONNAME2,
                 Distance_m = as.numeric(costDistance(trans, origin, goal)))
      } else {
        paths <-
          bind_rows(paths,
                    shortestPath(trans, origin, goal, output = "SpatialLines") %>% 
                      st_as_sf() %>% 
                      mutate(STATIONNAME1 = origin$STATIONNAME1,
                             STATIONNAME2 = origin$STATIONNAME2,
                             Distance_m = as.numeric(costDistance(trans, origin, goal))))
      }
      setTxtProgressBar(pb, i)
    }
  }, error=function(e){message("Error in calculating least cost path trajectories")})
  
  traj_lookup <- 
    traj_df %>% 
    dplyr::select(-c(Start.Longitude, Start.Latitude, 
                     End.Longitude, End.Latitude)) %>% 
    left_join(paths, by = c("STATIONNAME1", "STATIONNAME2")) %>% 
    st_sf()
  
  
  ## return list output with step distances and spatial trajectory file
  nonresidence <-
    residence$nonresidences %>% 
    dplyr::select(-c(DISTANCE, ROM, 
                     Start.Longitude, Start.Latitude, 
                     End.Longitude, End.Latitude)) %>% 
    left_join(traj_lookup, by = c("STATIONNAME1", "STATIONNAME2")) %>% 
    st_sf() %>% 
    mutate(ROM = Distance_m/DURATION)
  
  return(nonresidence)
}
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
traj_to_points <- function(traj, samp_interval = 30){
  
  traj <-
    traj %>% 
    mutate(num_pts = ceiling(traj$DURATION/(60*samp_interval)))
  
  traj_df <-
    traj %>% 
    as_tibble() %>% 
    dplyr::select(-geometry)
  
  pts <- 
    traj %>% 
    st_line_sample(n = traj$num_pts, type = "regular") %>% 
    st_sf() %>% 
    bind_cols(traj_df)
  
  return(pts)
  
}
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
## UD estimation
ud_est <- function(ATTdata, traj, residence, cost.ras, ll_epsg = 4326, timestep = 60,
                   utm_epsg, UDextent = 1, UDgrid = 50, h = 200, stepwise = FALSE){
  
  ll <- CRS(paste0("+init=epsg:", ll_epsg))
  utm <- CRS(paste0("+init=epsg:", utm_epsg))
  
  spdata_utm <-
    residence$residenceslog %>% 
    st_as_sf(coords = c("Station.Longitude", "Station.Latitude"), crs = ll_epsg) %>%
    st_transform(crs = utm_epsg) %>% 
    as_Spatial()
  
  UDwin <- owin(xrange = (extent(spdata_utm) + (diff(extent(spdata_utm)[1:2])*UDextent))[1:2],
                yrange = (extent(spdata_utm) + (diff(extent(spdata_utm)[3:4])*UDextent))[3:4])
  
  dimyx <- c(diff(UDwin$yrange)/UDgrid,
             diff(UDwin$xrange)/UDgrid)
  
  ## UD from COA data
  tryCatch({
    message("- Estimating Kernel Density from COA positions")
    
    coa.sf <-
      COA(ATTdata, timestep = timestep) %>%
      st_as_sf(coords=c("Longitude.coa","Latitude.coa"), crs = 4326, remove = F) %>%
      st_transform(utm_epsg)
    
    coa.snap <-
      coa.sf %>%
      as_Spatial() %>%
      maptools::snapPointsToLines(., as_Spatial(traj))
    
    coa.data <-coa.snap %>% as_tibble()
    
    pt.ppp <- ppp(x = coa.data$X, y = coa.data$Y, window = UDwin)
    pts.UD <- raster(density(pt.ppp, sigma = h, method = "C", dimyx = dimyx))
    values(pts.UD) <- abs(values(pts.UD)/max(values(pts.UD))-1) * 100
    projection(pts.UD) <- utm
    }, error=function(e){message("Error in estimating UD from COA positions")})
  
  ## UD from trajectory
  tryCatch({
    
    if(stepwise %in% TRUE){
      message("- Estimating stepwise Kernel Density from least cost path trajectory")
      ## Estimate sigma values for each movement step based on transit time
      sigs <- (abs(traj$DURATION/max(traj$DURATION)) * h) + h
      
      ## Estimate UD for each step
      for(l in 1:nrow(traj)){
        if(l %in% 1){
          t_stack <- stack()
          pb <- txtProgressBar(max = nrow(traj), style = 3)
        }
        step.traj <-
          traj %>%
          slice(l) %>%
          st_combine()
        step.psp <- as.psp(step.traj, window = UDwin)
        step.UD <- raster(density(step.psp, sigma = sigs[l], method = "C", dimyx = dimyx))
        projection(step.UD) <- utm
        t_stack <- addLayer(t_stack, step.UD)
        setTxtProgressBar(pb, l)
      }

      traj.UD <- sum(t_stack)
      values(traj.UD) <- abs(values(traj.UD)/max(values(traj.UD)) - 1) * 100
      
    } else {
      message("- Estimating Kernel Density from least cost path trajectory")
      trajectory_utm <-
        traj %>%
        st_combine()
      
      traj.psp <- as.psp(trajectory_utm, window = UDwin)
      traj.UD <- raster(density(traj.psp, sigma = h, method = "C", dimyx = dimyx))
      values(traj.UD) <- abs(values(traj.UD)/max(values(traj.UD)) - 1) * 100
      projection(traj.UD) <- utm 
    }
  }, error=function(e){message("Error in estimating UD for trajectories")})
  
  # Merge trajectory UD with point UD and remove land areas
  tryCatch({
    UDras <- mean(traj.UD, pts.UD)
    mask.ras <- resample(cost.ras, UDras, method = "ngb")
    values(mask.ras)[values(mask.ras) > 1] <- NA
    UDras.mask <- raster::mask(UDras, mask = mask.ras)
  }, error=function(e){message("Error in estimating merged UD")})
  
  # calculate UD areas
  tryCatch({
    message("- Estimating UD area for 50% and 95% contours")
    UD50 <- UDras < 50; UD50[values(UD50) %in% 0] <- NA
    UD50.area <- gArea(rasterToPolygons(UD50, dissolve=T))
    
    UD95 <- UDras < 95; UD95[values(UD95) %in% 0] <- NA
    UD95.area <- gArea(rasterToPolygons(UD95, dissolve=T))
  }, error=function(e){message("Error in estimating UD areas")})
  
  return(list(raster_output = UDras.mask,
              UD50_area = UD50.area,
              UD95_area = UD95.area))
}

## ************************************************************************************************************************************* ##
## adehabitat BBKUD formulation
# bbud_est <- function(pts, cost.ras, h = 200, ll_epsg = 4326, timestep = 60,
#                      utm_epsg, UDextent = 1, UDgrid = 50){
#  
#   require(adehabitatHR) 
#   
#   ll <- CRS(paste0("+init=epsg:", ll_epsg))
#   utm <- CRS(paste0("+init=epsg:", utm_epsg))
#   
#   dat <-
#     pts %>% 
#     st_cast("POINT") %>% 
#     group_by(TRANSMITTERID, NONRESIDENCEEVENT) %>% 
#     mutate(tt = seq(STARTTIME[1] + 1, ENDTIME[1] - 1, l = num_pts[1]))
#     
#   grext <- extent(dat) + (diff(extent(dat)[1:2])*UDextent)
#   gr <- expand.grid(x = seq(grext[1], grext[2], by = UDgrid), 
#                     y = seq(grext[3], grext[4], by = UDgrid))
#   coordinates(gr) <- ~x + y
#   projection(gr) <- utm
#   gridded(gr) <- TRUE
#   
#   tf <- as.ltraj(xy = st_coordinates(dat), date = dat$tt,
#                  id = dat$TRANSMITTERID, typeII = TRUE, proj4string = utm)
#   s1f <- (liker(tf, rangesig1 = c(0, 500), sig2 = h, 
#                 byburst = FALSE, plotit = FALSE)[[1]]$sig1)/div
#   tryCatch({
#     kbfull <- adehabitatHR::kernelbb(tf, sig1 = s1f, sig2 = h, grid = gr, byburst = T)
#     UDras <- raster(adehabitatHR::getvolumeUD(kbfull))
#   }, error = function(e) {
#     message("\nError in estimating merged UD")})
#   
#   mask.ras <- resample(cost.ras, UDras, method = "ngb")
#   values(mask.ras)[values(mask.ras) > 1] <- NA
#   UDras.mask <- raster::mask(UDras, mask = mask.ras)
#   
#   return(UDras.mask)
# }
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
### Function to plot output on leaflet
plot.lcUD <- function(UDobj, ...){
  library(mapview, warn.conflicts = F, quietly = T)
  library(leaflet, warn.conflicts = F, quietly = T)
  library(sf, warn.conflicts = F, quietly = T)
  library(raster, warn.conflicts = F, quietly = T)
  
  spdat <-
    UDobj$tagdata %>%  
    group_by(transmitter_id, station_name, longitude, latitude) %>% 
    summarise(.groups = 'keep',
              `Number of Detections` = n()) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  traj <- 
    UDobj$lc.traj %>% 
    mutate(transit_time_min = DURATION/60)

  UDras <- projectRaster(UDobj$UD.raster, crs = CRS("+init=epsg:4326"))
  UDras[raster::values(UDras) > 99] <- NA
  
  mapviewOptions(fgb = F)
  
  m <-
    mapview(UDras, layer = "Utilisation Distribution", homebutton = F, na.color = "transparent") +
    mapview(spdat, alpha = 0, layer = "Detection data", col.regions = "red", cex = "Number of Detections", legend = F, homebutton = F) +
    mapview(traj, zcol = "transit_time_min", layer = "Transit times (min)", homebutton = F)
    

  mm <-
    m@map %>% 
    leaflet::addLayersControl(baseGroups = c("CartoDB.Positron", "CartoDB.DarkMatter", "OpenStreetMap",
                                             "Esri.WorldImagery", "OpenTopoMap"),
                              overlayGroups = c("Utilisation Distribution", 
                                             "Detection data", "Transit times (min)"), 
                              options = leaflet::layersControlOptions(collapsed = F, position = "topleft")) %>% 
    setView(lat = mean(extent(UDras)[3:4]), lng = mean(extent(UDras)[1:2]), zoom = 10)
  
  mapviewOptions(fgb = T)
  
  return(mm)
}
# ## ************************************************************************************************************************************* ##






