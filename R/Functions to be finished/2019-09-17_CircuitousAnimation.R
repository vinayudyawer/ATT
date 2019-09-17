## Animating movements of animals along a circuitous linear system
## Written by Vinay Udyawer Sept 16 2019

library(tidyverse)
library(lubridate)
library(VTrack)
library(sf)
library(mapview)
library(move)
library(moveVis)

data("PointsCircuitous_crocs")
array <- as_tibble(PointsCircuitous_crocs)

data("crocs")
dets <-
  crocs %>%
  as_tibble() %>%
  left_join(array, by=c("Receiver.S.N" = "LOCATION")) %>%
  mutate(Date.Time = ymd_hms(Date.Time)) %>%
  arrange(ID, Date.Time) %>%
  st_as_sf(coords=c("LONGITUDE", "LATITUDE"), crs=4326) %>%
  group_by(ID) %>%
  mutate(Distance.between.fixes.meters = st_distance(geometry[lag(row_number())], geometry, by_element=T)) %>%
  as_Spatial() %>% as_tibble() %>%
  rename(LONGITUDE = coords.x1,
         LATITUDE = coords.x2)

detections <-
  dets %>%
  filter(Distance.between.fixes.meters > 0) %>%
  dplyr::select(-c(RADIUS, LONGITUDE, LATITUDE))

trajfun <-
  function(detections, array){
    combdat <-
      detections %>%
      left_join(array, by=c("Receiver.S.N" = "LOCATION"))

    for(i in 1:nrow(combdat)){
      if(i %in% 1){out <- combdat %>% slice(1)}

      tmpdf <- combdat %>% slice(i:(i+1))
      if(n_distinct(tmpdf$Receiver.S.N) > 1){

        waypoints <-
          array[which(grepl(tmpdf$Receiver.S.N[1], array$LOCATION)): which(grepl(tmpdf$Receiver.S.N[2], array$LOCATION)),] %>%
          slice(-1, -nrow(.))

        tmpdf <-
          tmpdf %>%
          add_row(LONGITUDE = waypoints$LONGITUDE,
                  LATITUDE = waypoints$LATITUDE,
                  .after = 1) %>%
          mutate(ID = first(ID),
                 Transmitter.Name = first(Transmitter.Name),
                 Transmitter.S.N = first(Transmitter.S.N),
                 Date.Time = seq(as.POSIXct(first(Date.Time)), as.POSIXct(last(Date.Time)), len = nrow(.)))
        out <- bind_rows(out, tmpdf)

        } else {
          out <- bind_rows(out, tmpdf)
        }
    }
    return(out %>% slice(-1))
  }

traj <-
  detections %>%
  group_by(ID) %>%
  do(trajfun(., array)) %>%
  filter(!duplicated(Date.Time, ID))

traj %>%
  st_as_sf(coords=c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
  mapview(zcol = "ID", burst = T)


## Use moveVis to animate movements
movedat <- df2move(df = traj,
                  proj = CRS("+init=epsg:4326"),
                  x = "LONGITUDE", y = "LATITUDE",
                  time = "Date.Time", track_id = "ID")

plot(movedat)

# align move_data to a uniform time scale
m <- align_move(movedat)
lines(m)

# create spatial frames
frames <-
  frames_spatial(m, map_service = "carto", map_type = "light") %>%
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_scalebar() %>%
  add_progress()

frames[[100]]

## This may take a while..
animate_frames(frames, fps = 25, end_pause = 2, out_file = "test/test.mp4")






