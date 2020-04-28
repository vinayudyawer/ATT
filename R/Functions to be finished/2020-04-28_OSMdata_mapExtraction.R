

get_map <- function(statinfo, expand = 1){
  require(tidyverse)
  require(raster)
  require(sf)
  require(osmdata)

  statinfo_sf <-
    statinfo %>%
    st_as_sf(coords = c("station_longitude", "station_latitude"),
             crs = 4326)

  bb <- extent(statinfo_sf) + expand

  polycoast <-
    opq(bbox = bb[c(1, 3, 2, 4)]) %>%
    add_osm_feature(key = 'natural', value = 'coastline') %>%
    osmdata_sf()

  islands <-
    polycoast$osm_polygons %>%
    rowid_to_column() %>%
    mutate(area = st_area(geometry)) %>%
    dplyr::select(rowid, area)

  coastline <- polycoast$osm_lines %>% st_union %>% st_line_merge
  pol <- st_as_sfc(st_bbox(bb)) %>% st_set_crs(4326)
  coastpoly <- lwgeom::st_split(st_geometry(pol), st_geometry(coastline))

  pol_list <- list()
  for (n in 1:length(coastpoly[[1]])) {
    pol_list[[n]] <-
      st_cast(coastpoly[[1]][n][[1]], 'POLYGON') %>% st_geometry() %>% st_sf()
    st_crs(pol_list[[n]]) <- 4326
  }

  land <-
    st_as_sf(data.table::rbindlist(pol_list)) %>%
    rowid_to_column() %>%
    slice(-2) %>%
    mutate(area = st_area(geometry)) %>%
    dplyr::select(rowid, area)

  studypol1 <-
    st_as_sf(data.table::rbindlist(list(land, islands), use.names = T)) %>%
    st_combine() %>%
    st_union()

  ## Remove freshwater
  polywater <-
    opq(bbox = bb[c(1, 3, 2, 4)]) %>%
    add_osm_feature(key = 'natural', value = "water") %>%
    osmdata_sf()

  polyriverbank <-
    opq(bbox = bb[c(1, 3, 2, 4)]) %>%
    add_osm_feature(key = 'waterway', value = "riverbank") %>%
    osmdata_sf()

  polyriver <-
    opq(bbox = bb[c(1, 3, 2, 4)]) %>%
    add_osm_feature(key = 'waterway', value = "river") %>%
    osmdata_sf()

  # rivers <-
  #   rbind(polyriver$osm_multilines %>% st_cast("LINESTRING") %>% transmute(rowid = osm_id),
  #         polyriver$osm_lines %>% transmute(rowid = osm_id)) %>%
  #   st_buffer(dist = 0.00001) %>%
  #   mutate(area = st_area(.))

  riverpol <-
    rbind(
      polywater$osm_polygons %>% filter(natural %in% "water") %>% transmute(rowid = osm_id, area = st_area(.)),
      polywater$osm_multipolygons %>% transmute(rowid = osm_id, area = st_area(.)),
      polyriverbank$osm_polygons %>% filter(waterway %in% "riverbank") %>% transmute(rowid = osm_id, area = st_area(.)),
      polyriverbank$osm_multipolygons %>% transmute(rowid = osm_id, area = st_area(.))
    ) %>%
    as_Spatial() %>%
    st_as_sf() %>%
    st_combine() %>%
    st_union()

  studymap <-
    studypol1 %>%
    st_buffer(dist = 0) %>%
    st_difference(riverpol %>% st_buffer(dist = 0)) %>%
    as_Spatial() %>%
    st_as_sf

  return(studymap)
}











