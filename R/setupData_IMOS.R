#' Setup data prior to use with Animal Tracking Toolbox functions
#'
#'@description Sets up and stores all relevant information required for spatio-temporal analyses of passive telemetry data.
#'  This function produces an 'ATT' object which standardises field names (currently handles IMOS ATF export data structure).
#'
#' @param Tag.Detections data frame with tag detection data including coordinates. Currently only handles IMOS ATF data structure.
#' @param Tag.Metadata data frame with metadata for all tags represented in Tag.Detections. Currently only handles IMOS ATF metadata structure.
#' @param Station.Information data frame with information on receiver station including coordinates. Currently only handles IMOS ATF station information structure.
#' @param tz time zone of date time information in Tag.Detections, Tag.Metadata and Station.Information. If none provided defaults to "UTC"
#' @param crs geographic coordinate system for all Tag.Detections, Tag.Metadata and Station.Information (latitude/longitude). If none provided defaults to WGS84.
#'
#' @return Produces an 'ATT' object that is a list of tibbles containing Tag.Detections, Tag.Metadata and Station.Information.
#'   The 'ATT' object will have a geographic coordinate system associated with it for smoother functioning of subsequent functions.
#'
#' @seealso setup data can be used to estimate detection \code{\link{detectionSummary}}, dispersal \code{\link{dispersalSummary}}
#'   and Short-term center of activity \code{\link{COA}}.
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr select
#' @importFrom lubridate ymd_hms
#' @importFrom lubridate date
#' @importFrom sp CRS
#' @examples
#' ## Import example datasets
#' data(tagdata)
#' data(taginfo)
#' data(statinfo)
#'
#' ## Setup data
#' ATTdata<- setupData(Tag.Detections = tagdata, Tag.Metadata = taginfo, Station.Information = statinfo)
#' ATTdata
#'
#'
setupData_IMOS<-function(Tag.Detections, Tag.Metadata, Station.Information, tz="UTC", crs=NULL){

  object<-
    structure(
      list(
        Tag.Detections = as_tibble(Tag.Detections) %>%
          transmute(Date.Time = lubridate::ymd_hms(detection_timestamp, tz=tz),
                    Transmitter = transmitter_id,
                    Station.Name = station_name,
                    Receiver = receiver_name,
                    Latitude = latitude,
                    Longitude = longitude,
                    Sensor.Value = sensor_value,
                    Sensor.Unit = sensor_unit),

        Tag.Metadata = as_tibble(Tag.Metadata) %>%
          transmute(Tag.ID = tag_id,
                    Transmitter = transmitter_id,
                    Sci.Name = scientific_name,
                    Common.Name = common_name,
                    Tag.Project = tag_project_name,
                    Release.Latitude = release_latitude,
                    Release.Longitude = release_longitude,
                    Release.Date = lubridate::date(ymd_hms(ReleaseDate, tz=tz)),
                    Tag.Life = tag_expected_life_time_days,
                    Tag.Status = tag_status,
                    Sex = sex,
                    Bio = measurement),

        Station.Information = as_tibble(Station.Information) %>%
          transmute(Station.Name = station_name,
                    Receiver = receiver_name,
                    Installation = installation_name,
                    Receiver.Project = project_name,
                    Deployment.Date = lubridate::ymd_hms(deploymentdatetime_timestamp, tz=tz),
                    Recovery.Date = lubridate::ymd_hms(recoverydatetime_timestamp, tz=tz),
                    Station.Latitude = station_latitude,
                    Station.Longitude = station_longitude,
                    Receiver.Status = status)),

      class="ATT")

  if(inherits(crs, "CRS")){
    attr(object, "CRS")<-crs
    }else{
      message("Geographic projection for detection positions not recognised, reverting to WGS84 global coordinate reference system")
      attr(object, "CRS")<-CRS("+init=epsg:4326")
    }

  return(object)
}

