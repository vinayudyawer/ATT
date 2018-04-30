
COA<-function (ATTdata, id="Tag.ID", timestep=60, split=FALSE){

  suppressPackageStartupMessages({require(dplyr); require(lubridate)})

  if(!inherits(ATTdata, "ATT"))
    stop("Oops! Input data needs to be an 'ATT' object.
         \nSet up your data first using setupData() before running this operation")

  ## Combine Tag.Detection and Tag.Metadata into a combined tibble for processing
  data<- left_join(ATTdata$Tag.Detections, ATTdata$Tag.Metadata, by="Tag.ID") %>%
    mutate_at(id, factor)

  step_sec <- timestep * 60
  ex <- seq(from = trunc(min(data$Date.Time, na.rm = TRUE), "day"),
            to = trunc(max(data$Date.Time, na.rm = TRUE), "day") + 86400,
            by = step_sec)
  data$TimeStep.coa <- cut(data$Date.Time, breaks = ex)

  cenac <-
    data %>%
    group_by_at(vars(id, TimeStep.coa)) %>%
    summarize(Latitude.coa = mean(Latitude, na.rm = TRUE),
              Longitude.coa = mean(Longitude, na.rm = TRUE),
              Sensor.Value.coa = mean(Sensor.Value),
              Sensor.Unit = first(Sensor.Unit),
              Number.of.Stations = n_distinct(Station.Name),
              Number.of.Detections = n(),
              Sci.Name = first(Sci.Name),
              Common.Name = first(Common.Name),
              Tag.Project = first(Tag.Project),
              Release.Latitude = first(Release.Latitude),
              Release.Longitude = first(Release.Longitude),
              Release.Date = first(Release.Date),
              Tag.Life = first(Tag.Life),
              Tag.Status = first(Tag.Status),
              Sex = first(Sex),
              Bio = first(Bio)) %>%
    mutate(TimeStep.coa = lubridate::ymd_hms(TimeStep.coa))

  if(length(group_size(cenac)) > 1 & split == TRUE){
    cenac <- split(cenac, cenac[,id])
    attr(cenac, "class")<-c("list","COA", "ATT")
    attr(cenac, "CRS")<-attr(ATTdata, "CRS")
  }else{
    attr(cenac, "class")<-c("tbl_df","tbl","data.frame","COA", "ATT")
    attr(cenac, "CRS")<-attr(ATTdata, "CRS")
  }
  return(cenac)
}


