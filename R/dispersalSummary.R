
dispersalSummary<-function(ATTdata){

  suppressPackageStartupMessages({require(dplyr); require(sp); require(raster); require(maptools)})

  if(!inherits(ATTdata, "ATT"))
    stop("Oops! Input data needs to be an 'ATT' object.
         \nSet up your data first using setupData() before running this operation")

  ## Combine Tag.Detection and Tag.Metadata into a combined tibble for processing
  data<- left_join(ATTdata$Tag.Detections, ATTdata$Tag.Metadata, by="Tag.ID")

  ## CRS for geographic coordinates
  ll<-CRS("+init=epsg:4326") ## Geographic projection for lat/long data

  dispfun<-function(dat){
    if(!is.na(dat$Release.Latitude[1])){
      pts<-data.frame(lat=dat$Latitude, lon=dat$Longitude); coordinates(pts)<-~lon+lat; projection(pts)<-ll
      pt<-data.frame(lat=dat$Release.Latitude[1], lon=dat$Release.Longitude[1]); coordinates(pt)<-~lon+lat; projection(pt)<-ll

      disp<-
        dat%>%
        ### Straight line distance between release location and each detection
        mutate(Release.Dispersal = spDistsN1(pts=pts, pt=pt, longlat=TRUE))
      disp$Release.Bearing<- sapply(1:nrow(disp),
                                   function(x)
                                     gzAzimuth(from=matrix(c(disp$Release.Longitude[x],disp$Release.Latitude[x]),1,2),
                                               to=matrix(c(disp$Longitude[x],disp$Latitude[x]),1,2)))
      disp$Release.Bearing<- ifelse(!is.na(disp$Release.Bearing) & disp$Release.Bearing<0, disp$Release.Bearing+360, disp$Release.Bearing)

      ### Straight line distance between consecutive detections
      disp$Consecutive.Dispersal<- c(spDistsN1(matrix(c(disp$Release.Longitude[1],disp$Release.Latitude[1]),1,2),
                                               matrix(c(disp$Longitude[1],disp$Latitude[1]),1,2), longlat=TRUE),
                                     sapply(2:nrow(disp),
                                            function(x)
                                              spDistsN1(matrix(c(disp$Longitude[x-1],disp$Latitude[x-1]),1,2),
                                                        matrix(c(disp$Longitude[x],disp$Latitude[x]),1,2), longlat=TRUE)))
      disp$Consecutive.Bearing<- c(gzAzimuth(from=matrix(c(disp$Release.Longitude[1],disp$Release.Latitude[1]),1,2),
                                             to=matrix(c(disp$Longitude[1],disp$Latitude[1]),1,2)),
                                   sapply(2:nrow(disp),
                                          function(x)
                                            gzAzimuth(from=matrix(c(disp$Longitude[x-1],disp$Latitude[x-1]),1,2),
                                                      to=matrix(c(disp$Longitude[x],disp$Latitude[x]),1,2))))
      disp$Consecutive.Bearing<- ifelse(!is.na(disp$Consecutive.Bearing) & disp$Consecutive.Bearing<0, disp$Consecutive.Bearing+360, disp$Consecutive.Bearing)
    }
    else{
      pts<-data.frame(lat=dat$Latitude, lon=dat$Longitude); coordinates(pts)<-~lon+lat; projection(pts)<-ll

      disp<-
        dat%>%
        ### Straight line distance between release location and each detection
        mutate(Release.Dispersal = NA,
               Release.Bearing = NA)

      ### Straight line distance between consecutive detections
      disp$Consecutive.Dispersal<- c(NA, sapply(2:nrow(disp),
                                                function(x)
                                                  spDistsN1(matrix(c(disp$Longitude[x-1],disp$Latitude[x-1]),1,2),
                                                            matrix(c(disp$Longitude[x],disp$Latitude[x]),1,2), longlat=TRUE)))
      disp$Consecutive.Bearing<- c(NA, sapply(2:nrow(disp),
                                              function(x)
                                                gzAzimuth(from=matrix(c(disp$Longitude[x-1],disp$Latitude[x-1]),1,2),
                                                          to=matrix(c(disp$Longitude[x],disp$Latitude[x]),1,2))))
      disp$Consecutive.Bearing<- ifelse(!is.na(disp$Consecutive.Bearing) & disp$Consecutive.Bearing<0, disp$Consecutive.Bearing+360, disp$Consecutive.Bearing)
    }
    return(disp)
  }


  disptab<-
    data %>%
    group_by(Tag.ID) %>%
    do(dispfun(.)) %>%
    ungroup()

  return(disptab)
}
