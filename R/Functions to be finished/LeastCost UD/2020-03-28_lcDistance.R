#### Least cost path calculation using passive telemetry datasets
## Written by VU 2018-02-12
## Last updated: 2020-03-28 (VU)
##
## Function to calculate least cost paths between detections from an ATTdata object (from VTrack)
## around landmasses. Uses 'osmdata' to construct cost raster if none is provided

## ************************************************************************************************************************************* ##
lcDistance <- function(ATTdata, cost = NULL, trans = NULL, utm_epsg, ll_epsg = 4326, timestep = 60, h = 100, UDextent = 1, UDgrid = 100, cost.res = 100, directions = 16, ...){
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
  sapply(c("lubridate","sf","dplyr","raster","gdistance","spatstat","maptools","rgeos","VTrack", "data.table", "lwgeom"), require, character.only=TRUE, warn.conflicts = F, quietly = T)
  ll <- CRS(paste0("+init=epsg:", ll_epsg))
  utm <- CRS(paste0("+init=epsg:", utm_epsg))

  combdata <-
    ATTdata$Tag.Detections %>%
    left_join(ATTdata$Station.Information %>%
                dplyr::select(Receiver, Station.Name, Station.Latitude, Station.Longitude),
              by=c("Station.Name", "Receiver"))

  statinfo <-
    ATTdata$Station.Information %>%
    st_as_sf(coords=c("Station.Longitude", "Station.Latitude"), crs = ll_epsg)

  ## Setup input tagdata; convert to spatial object and transform to projection
  spdata_ll <-
    combdata %>%
    filter(!is.na(Station.Longitude)) %>%
    mutate(step = 0:(nrow(.)-1)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = ll_epsg)

  spdata_utm <-
    spdata_ll %>%
    st_transform(crs = utm_epsg)

  ## Extract landmass to calculate transition and cost raster if not provided
  if(is.null(trans)){
    if(!is.null(cost)){
      message("\n- No transition layer provided. Accessing cost layer provided")
      cost.in_utm <- raster::projectRaster(cost, crs=utm, method = "ngb")
      cost.ras <- resample(cost.in_utm, raster(extent(cost.in_utm), res = cost.res), method = "ngb")
      projection(cost.ras) <- utm

      ## Produce transition matrices, and correct for distortion
      tryCatch({
        message("\n- Calculating transition matrices for least cost path estimation")
        trCost <- transition(1/cost.ras, mean, directions = directions)
        trCost <- geoCorrection(trCost, type = "c")
      }, error=function(e){message("\nError in calculating Transition layer")})

    } else {
      message("\n- No cost or transition layer provided. Downloading coastline data from Open Street Map server")
      require("osmdata")
      tryCatch({

        bb <- extent(statinfo) + 0.1

        polydat <-
          opq(bbox = bb[c(1,3,2,4)]) %>%
          add_osm_feature(key = 'natural', value = 'coastline') %>%
          osmdata_sf

        islands <-
          polydat$osm_polygons %>%
          rowid_to_column() %>%
          mutate(area = st_area(geometry)) %>%
          dplyr::select(rowid, area)

        coastline <- polydat$osm_lines %>% st_union %>% st_line_merge
        pol <- st_as_sfc(st_bbox(coastline))
        coastpoly <- lwgeom::st_split(st_geometry(pol), st_geometry(coastline))

        pol_list <- list()
        for(n in 1:length(coastpoly[[1]])){
          pol_list[[n]] <- st_cast(coastpoly[[1]][n][[1]], 'POLYGON') %>% st_geometry() %>% st_sf()
          st_crs(pol_list[[n]]) <- 4326
        }

        land <-
          st_as_sf(data.table::rbindlist(pol_list)) %>%
          rowid_to_column() %>%
          mutate(area = st_area(geometry)) %>%
          filter(area %in% max(area))

        studypol <- st_as_sf(data.table::rbindlist(list(land, islands), use.names = T))

        poly <-
          studypol %>%
          st_transform(utm_epsg) %>%
          st_crop(extent(statinfo %>% st_transform(utm_epsg)) + 5000)

        cost.ras<-rasterize(poly, raster(extent(statinfo %>% st_transform(utm_epsg)) + 5000, res = cost.res), 1000)
        cost.ras[is.na(values(cost.ras))] <- 1
        projection(cost.ras) <- utm

        ## Produce transition matrices, and correct for distortion
        tryCatch({
          message("\n- Calculating transition matrices for least cost path estimation")
          trCost <- transition(1/cost.ras, mean, directions = directions)
          trCost <- geoCorrection(trCost, type = "c")
        }, error=function(e){message("\nError in calculating Transition layer")})

      }, error=function(e){message("\nError: no land in sight!\nConsider adding your own cost layer")})
      }
    } else {
      message("\n- Accessing provided transition layer for least cost path estimation")
      trCost <- trans
    }

  ## Construct shortest path trajectories between sequence of detection steps
  outdat <-
    spdata_ll %>%
    as_Spatial %>%
    as_tibble %>%
    dplyr::rename(Longitude = coords.x1, Latitude = coords.x2) %>%
    mutate(Transit.time_sec = as.numeric(difftime(Date.Time, lag(Date.Time), units = "sec")),
           Distance_m = NA) %>%
    dplyr::select(-step)

  tryCatch({
    message("- Constructing least cost path trajectories between consecutive detections")
    traj <- list()
    for(i in 1:max(spdata_utm$step)){
      if(i %in% 1){pb <- txtProgressBar(min=1, max=max(spdata_utm$step), style=3)}

      stepdat <- spdata_utm %>% filter(step %in% c(i-1, i))
      origin <- stepdat %>% slice(1) %>% as_Spatial
      goal <- stepdat %>% slice(n()) %>% as_Spatial

      if (nrow(distinct(stepdat)) > 1){
        traj[[i]] <- shortestPath(trCost, origin, goal, output = "SpatialLines")
        outdat$Distance_m[i+1] <- as.numeric(costDistance(trCost, origin, goal))
      } else {
        traj[[i]]<- NULL
        outdat$Distance_m[i+1]<-0
      }
      setTxtProgressBar(pb, i)
    }
    
    traj_clean <- Filter(Negate(is.null), traj)
    trajectory_utm <- do.call(rbind, traj_clean)
    trajectory_ll <- 
      st_as_sf(trajectory_utm) %>% 
      st_transform(ll_epsg) %>% 
      mutate(Transit.time_min = outdat %>% filter(Distance_m > 0) %>% .$Distance_m/60)
  }, error=function(e){message("\nError in calculating least cost path trajectories")})

  ## Calculate Kernel Density from trajectory
  UDwin <- owin(xrange = (extent(spdata_utm) + (diff(extent(spdata_utm)[1:2])*UDextent))[1:2],
                yrange = (extent(spdata_utm) + (diff(extent(spdata_utm)[3:4])*UDextent))[3:4])
  dimyx <- c(diff(UDwin$yrange)/UDgrid,
             diff(UDwin$xrange)/UDgrid)

  # point UD
  tryCatch({
    message("\n- Estimating Kernel Density from COA positions")
    coa.sf <-
      COA(ATTdata, timestep = timestep) %>%
      st_as_sf(coords=c("Longitude.coa","Latitude.coa"), crs = 4326, remove = F) %>%
      st_transform(utm_epsg)
    suppressWarnings(coa.snap <- 
                       coa.sf %>% 
                       as_Spatial() %>% 
                       snapPointsToLines(., trajectory_utm))
    coa.data <-coa.snap %>% as_tibble()
    suppressWarnings(pt.ppp <- ppp(x = coa.data$X, y = coa.data$Y, window = UDwin))
    suppressWarnings(UDras.pt <- raster(density(pt.ppp, sigma = h, method = "C", dimyx = dimyx)))
    values(UDras.pt) <- abs(values(UDras.pt)/max(values(UDras.pt)) - 1) * 100
    projection(UDras.pt) <- utm
  }, error=function(e){message("Error in estimating UD from COA positions")})

  # trajectory UD
  tryCatch({
    message("- Estimating Kernel Density from least cost path trajectory")
    traj.psp <- as.psp(trajectory_utm, window = UDwin)
    ## Estimate sigma values for each movement step based on transit time
    # sigs <- (abs(trajectory_ll$Transit.time_min/max(trajectory_ll$Transit.time_min)) * h)
    ## Estimate UD for each step
    # for(l in 1:length(trajectory_utm)){
    #   if(l %in% 1){
    #     t_stack <- stack()
    #     pb <- txtProgressBar(max = length(trajectory_utm), style = 3)
    #   }
    #   step.psp <- as.psp(trajectory_utm[l], window = UDwin)
    #   step.UD <- raster(density(step.psp, sigma = sigs[l], method = "C", dimyx = dimyx))
    #   values(step.UD) <- abs(values(step.UD)/max(values(step.UD)) - 1) * 100
    #   projection(step.UD) <- utm
    #   t_stack <- addLayer(t_stack, step.UD)  
    #   setTxtProgressBar(pb, l)
    # }
    UDras.traj <- raster(density(traj.psp, sigma = h, method = "C", dimyx = dimyx))
    values(UDras.traj) <- abs(values(UDras.traj)/max(values(UDras.traj)) - 1) * 100
    projection(UDras.traj) <- utm
    # UDras.traj <- mean(t_stack)
  }, error=function(e){message("\nError in estimating UD for trajectories")})

  # Merge trajectory UD with point UD and remove land areas
  tryCatch({
    UDras <- mean(UDras.traj, UDras.pt)
    mask.ras <- resample(cost.ras, UDras, method = "ngb")
    values(mask.ras)[values(mask.ras) > 1] <- NA
    UDras.m_utm <- raster::mask(UDras, mask = mask.ras)
    UDras.m_ll <- projectRaster(UDras.m_utm, crs = ll)
  }, error=function(e){message("Error in estimating UD")})

  # calculate UD areas
  tryCatch({
    message("\n- Estimating UD area for 50% and 95% contours")
    UD50 <- UDras.m_utm < 50; UD50[values(UD50) %in% 0] <- NA
    UD50.area <- gArea(rasterToPolygons(UD50, dissolve=T))
    
    UD95 <- UDras.m_utm < 95; UD95[values(UD95) %in% 0] <- NA
    UD95.area <- gArea(rasterToPolygons(UD95, dissolve=T))
    
    }, error=function(e){message("Error in estimating UD areas")})

  ## return list output with step distances and spatial trajectory file
  out<-list(tagdata = outdat,
            coa.data = coa.sf,
            kernel.areas = data.frame(UD50_m2 = UD50.area, UD95_m2 = UD95.area),
            lc.traj = trajectory_ll,
            UD.raster = UDras.m_ll,
            cost.raster = projectRaster(cost.ras, crs = ll, method = "ngb"))

  return(out)
}
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
### Function to plot output on leaflet
plot.lcUD <- function(UDobj, ...){
  library(mapview, warn.conflicts = F, quietly = T)
  library(sf, warn.conflicts = F, quietly = T)
  library(raster, warn.conflicts = F, quietly = T)
  spdat <-
    UDobj$tagdata %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    group_by(Station.Name) %>%
    summarise(`Number of Detections` = n())

  UDras <- UDobj$UD.raster
  UDras[raster::values(UDras) > 99] <- NA
  # coadat <-
  #   UDobj$coa.data %>%
  #   st_as_sf(coords = c("Longitude.coa", "Latitude.coa"), crs = 4326)

  m <-
    mapview(UDras, layer = "Utilisation Distribution", homebutton = F, na.color = "transparent", ...) +
    mapview(spdat, alpha = 0, layer = "Detection data", col.regions = "red", cex = "Number of Detections", legend = F, homebutton = F) +
    mapview(UDobj$lc.traj, zcol = "Transit.time_min", layer = "Transit times (min)", homebutton = F)

  return(m)
}
## ************************************************************************************************************************************* ##


## ************************************************************************************************************************************* ##
## Function to generate least cost trajectories
lcTraj <- function(ATTdata, cost = NULL, utm_epsg, ll_epsg = 4326, cost.res = 50, directions = 16, ...){

  ## load required libraries and set up CRSs
  sapply(c("lubridate","sf","dplyr","raster","gdistance","spatstat","maptools","rgeos","VTrack", "data.table", "lwgeom"), require,
         character.only=TRUE, warn.conflicts = F, quietly = T)
  ll <- CRS(paste0("+init=epsg:", ll_epsg))
  utm <- CRS(paste0("+init=epsg:", utm_epsg))

  combdata <-
    ATTdata$Tag.Detections %>%
    dplyr::select(-c(Longitude, Latitude)) %>%
    left_join(ATTdata$Station.Information %>%
                dplyr::select(Receiver, Station.Name, Station.Latitude, Station.Longitude),
              by=c("Station.Name", "Receiver")) %>%
    rename(Longitude = Station.Longitude,
           Latitude = Station.Latitude)

  statinfo <-
    ATTdata$Station.Information %>%
    st_as_sf(coords=c("Station.Longitude", "Station.Latitude"), crs = ll_epsg)

  ## Setup input tagdata; convert to spatial object and transform to projection
  spdata_ll <-
    combdata %>%
    mutate(step = 0:(nrow(.)-1)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = ll_epsg)

  spdata_utm <-
    spdata_ll %>%
    st_transform(crs = utm_epsg)

  ## Extract landmass to calculate transition and cost raster if not provided
  if(is.null(trans)){
    if(!is.null(cost)){
      message("\n- No transition layer provided. Accessing cost layer provided")
      cost.in_utm <- raster::projectRaster(cost, crs=utm, method = "ngb")
      cost.ras <- resample(cost.in_utm, raster(extent(cost.in_utm), res = cost.res), method = "ngb")
      projection(cost.ras) <- utm
      
      ## Produce transition matrices, and correct for distortion
      tryCatch({
        message("\n- Calculating transition matrices for least cost path estimation")
        trCost <- transition(1/cost.ras, mean, directions = directions)
        trCost <- geoCorrection(trCost, type = "c")
      }, error=function(e){message("\nError in calculating Transition layer")})
      
    } else {
      message("\n- No cost or transition layer provided. Downloading coastline data from Open Street Map server")
      require("osmdata")
      tryCatch({
        
        bb <- extent(statinfo) + 0.1
        
        polydat <-
          opq(bbox = bb[c(1,3,2,4)]) %>%
          add_osm_feature(key = 'natural', value = 'coastline') %>%
          osmdata_sf
        
        islands <-
          polydat$osm_polygons %>%
          rowid_to_column() %>%
          mutate(area = st_area(geometry)) %>%
          dplyr::select(rowid, area)
        
        coastline <- polydat$osm_lines %>% st_union %>% st_line_merge
        pol <- st_as_sfc(st_bbox(coastline))
        coastpoly <- lwgeom::st_split(st_geometry(pol), st_geometry(coastline))
        
        pol_list <- list()
        for(n in 1:length(coastpoly[[1]])){
          pol_list[[n]] <- st_cast(coastpoly[[1]][n][[1]], 'POLYGON') %>% st_geometry() %>% st_sf()
          st_crs(pol_list[[n]]) <- 4326
        }
        
        land <-
          st_as_sf(data.table::rbindlist(pol_list)) %>%
          rowid_to_column() %>%
          mutate(area = st_area(geometry)) %>%
          filter(area %in% max(area))
        
        studypol <- st_as_sf(data.table::rbindlist(list(land, islands), use.names = T))
        
        poly <-
          studypol %>%
          st_transform(utm_epsg) %>%
          st_crop(extent(statinfo %>% st_transform(utm_epsg)) + 5000)
        
        cost.ras<-rasterize(poly, raster(extent(statinfo %>% st_transform(utm_epsg)) + 5000, res = cost.res), 1000)
        cost.ras[is.na(values(cost.ras))] <- 1
        projection(cost.ras) <- utm
        
        ## Produce transition matrices, and correct for distortion
        tryCatch({
          message("\n- Calculating transition matrices for least cost path estimation")
          trCost <- transition(1/cost.ras, mean, directions = directions)
          trCost <- geoCorrection(trCost, type = "c")
        }, error=function(e){message("\nError in calculating Transition layer")})
        
      }, error=function(e){message("\nError: no land in sight!\nConsider adding your own cost layer")})
    }
  } else {
    message("\n- Accessing provided transition layer for least cost path estimation")
    trCost <- trans
  }

  ## Construct shortest path trajectories between sequence of detection steps
  outdat <-
    spdata_ll %>%
    as_Spatial %>%
    as_tibble %>%
    rename(Longitude = coords.x1, Latitude = coords.x2) %>%
    mutate(Distance_m = NA) %>%
    dplyr::select(-step)

  tryCatch({
    message("- Constructing least cost path trajectories between consecutive detections")
    traj<-list()
    for(i in 1:max(spdata_utm$step)){
      if(i %in% 1){pb <- txtProgressBar(min=1, max=max(spdata_utm$step), style=3)}

      stepdat <- spdata_utm %>% filter(step %in% c(i-1, i))
      origin <- stepdat %>% slice(1) %>% as_Spatial
      goal <- stepdat %>% slice(n()) %>% as_Spatial

      if (nrow(distinct(stepdat)) > 1){
        traj[[i]]<-shortestPath(trCost, origin, goal, output = "SpatialLines")
        outdat$Distance_m[i+1]<-as.numeric(costDistance(trCost, origin, goal))
      } else {
        traj[[i]]<- NULL
        outdat$Distance_m[i+1]<-0
      }
      setTxtProgressBar(pb, i)
    }

    traj_clean <- Filter(Negate(is.null), traj)
    trajectory_utm <- do.call(rbind, traj_clean)
    trajectory_ll <- 
      st_as_sf(trajectory_utm) %>% 
      st_transform(ll_epsg) %>% 
      mutate(Transit.time_min = outdat %>% filter(Distance_m > 0) %>% .$Distance_m/60)
  }, error=function(e){message("\nError in calculating least cost path trajectories")})

  ## return list output with step distances and spatial trajectory file
  out<-list(tagdata = outdat,
            lc.traj = trajectory_ll,
            cost.raster = projectRaster(cost.ras, crs = ll, method = "ngb"))


  return(out)
}
## ************************************************************************************************************************************* ##



