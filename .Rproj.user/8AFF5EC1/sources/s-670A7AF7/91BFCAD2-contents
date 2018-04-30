### Animal Tracking Toolbox
###
### script written by Vinay Udyawer and the IMOS Task Team 
### Function uses raw detection data from passive telemetry or node-based datasets to summarise movement and acitivity space metrics
### Movement and space use metrics are calculated over full taglife and user defined temporal subsets


ATT<-function(tagdata, taginfo, IMOSdata=FALSE, sig2=200, timestep=60, ext=2, grid=200, sub="%Y-%m", cumulative=FALSE, plotfull=FALSE, plotsub=FALSE, storepoly=FALSE, QCval=2, div=4){
  ####################################################################################################################################################
  ### tagdata:    detection data for each indivdual (sourced from IMOS ATF or from VUE)
  ### taginfo:    Metadata information about each tag (e.g. sex, size, transmitter life, etc.)
  ### IMOSdata:   Flag to identify if data is from IMOS/AODN database; if FALSE then standard VUE output accepted
  ### sig2:       Smoothing factor used for Brownian Bridge KUD (in m), related to imprecission of relocations [default= 200 m for acoustic recievers]
  ### timestep:   Timestep used for center of activity subsetting (in minutes) [default= 60 min]
  ### ext:        Extent metric used to define extent of grid for brownian bridge kud estimation [default= 2]
  ### grid:       Grid resolution width for brownian bridge kud estimation [default= 200]
  ### sub:        Level of subsampling ** smaller the subsample = less data for KUD calculations
  ###               yearmon   = "%Y-%m" [default]
  ###               yearweek  = "%Y-%W"
  ###               yearly    = "%Y"
  ###               monthly   = "%m"
  ###               weekly    = "%W"
  ### cumulative: Calculate cumulative MCP and KUD area values **Time Intensive** [default= FALSE]
  ### plotfull:   Plot full bbKUD estimates with centre of activity points [default= FALSE]
  ### plotsub:    Plot subsetted bbKUD estimates with centre of activity points [default= FALSE]
  ### storepoly:  Store polygons of mcp, 50% and 95% full bbKUD as shapefiles in resulting list [default= FALSE]
  ### Qval:       Maximum Quality control flag to use when using IMOS ATF data 1:valid detection, 2:probably valid detection, 3:probably bad detection, 4:bad detection) [default= 2]
  ### div:        Divisor for sig1 smoothing factor [default = 4]          
  ###
  ####################################################################################################################################################
  ##### Setup libraries
  sapply(c("maptools", "maps", "adehabitatHR", "sp", "raster", "rasterVis", "plyr", "lubridate", "rgeos"), require, character.only=TRUE)
  
  
  ##### Setup input data based on source (IMOS or VEMCO)
  if(IMOSdata){data<-merge(subset(tagdata, Detection_QC<=QCval),taginfo, by="transmitter_id")}else{
    data<-merge(data.frame(tag_id=tagdata$Transmitter.Serial, transmitter_id=tagdata$Transmitter, 
                           installation_name=NA,
                           station_name=tagdata$Station.Name,
                           receiver_name=tagdata$Receiver,
                           detection_timestamp=tagdata[,grep("*Date*",colnames(tagdata))],
                           longitude=tagdata$Longitude,
                           latitude=tagdata$Latitude,
                           sensor_value=tagdata$Sensor.Value,
                           sensor_unit=tagdata$Sensor.Unit), taginfo, by="transmitter_id")
  }
  
  if(nrow(data)>1){
    message('Animal Tracking Toolbox: Calculating standardised metrics for tag ', as.character(data$transmitter_id[1]) )
    ##### Setup data
    ll<-CRS("+proj=longlat +datum=WGS84"); utm<-CRS("+init=epsg:3577")
    data$dt<-ymd_hms(data$detection_timestamp, tz="GMT")
    data$yearmon<-as.factor(format(data$dt, sub))
    sex<-data$sex[1]
    taglife<-data$tag_expected_life_time_days[1]
    bio<-data$measurement[1]
    
    ##### List preprocessing
    if(nrow(data)>1){
      #### Summarising data for overall ('full') and subsetted ('res')
      full<-data.frame(tag_id=data$tag_id[1], transmitter_id=data$transmitter_id[1],
                       species=data$scientific_name[1], 
                       sex=sex, bio=bio,
                       num_det=sum(table(as.Date(data$dt))[table(as.Date(data$dt))>2]),
                       days_det=length(table(as.Date(data$dt))[table(as.Date(data$dt))>2]), ## number of days detected, with >2 detections per day
                       num_stat=length(unique(data$station_name)))
      res<-data.frame(tag_id=data$tag_id[1], transmitter_id=data$transmitter_id[1],
                      species=data$scientific_name[1], 
                      sex=sex, bio=bio,
                      yearmon=data.frame(table(data$yearmon))[,1], 
                      num_det=data.frame(table(data$yearmon))[,2],
                      days_det=aggregate(as.Date(dt)~yearmon, data, function(x) length(table(x)[table(x)>2]), na.action=NULL)[,2], ## number of days detected, with >2 detections per day
                      num_stat=aggregate(station_name~yearmon, data, function(x) length(unique(x)), na.action=NULL)[,2],
                      num_new_stat=NA)
      ### Number of new stations detected on between timesteps
      if(nrow(res)>1){res$num_new_stat<-c(NA,sapply(2:nrow(res), function(x) num_new_stat=length(setdiff(unique(data[data$yearmon%in%res$yearmon[x],"station_name"]),unique(data[data$yearmon%in%res$yearmon[x-1],"station_name"])))))}else{res$num_new_stat<-NA}
      
      ### Calculate and add Detection Index (Residency Index: Number of days detected on IMOS array/ number of days in month)
      ### Need information on tag life to refine the index for first and last month of results
      if(!is.na(as.Date(data$ReleaseDate)[1])){
        st<-as.Date(data$ReleaseDate)[1] ## start of tag transmission (presuming release date)
      }else{
        st<-as.Date(min(data$dt, na.rm=TRUE)) ## use first date of detection as start date if release date not provided
      }
      if(!is.na(taglife)){
        et<-st+taglife ## maximum tag life from metadata to find date when battery life ends
        full$DI<-full$days_det/taglife
      }else{
        et<-as.Date(max(data$dt, na.rm=TRUE)) ## if not recorded in metadata, last date of detection used as end date
        full$DI<-full$days_det/as.numeric(difftime(et, st,"days")) ## DI when tag life not recorded (using last day detected as end date)
      }
      dal<-data.frame(yearmon=res$yearmon, dal=days_in_month(as.POSIXct(ymd(paste(res$yearmon,01,sep="-")), tz="GMT")))
      dal[dal$yearmon%in%format(st, "%Y-%m"),"dal"]<-dal[dal$yearmon%in%format(st, "%Y-%m"),"dal"]-as.integer(format(st, "%d"))
      dal[dal$yearmon%in%format(et, "%Y-%m"),"dal"]<-as.integer(format(et, "%d"))
      dal$DI<-res$days_det/dal$dal
      res<-merge(res, dal[,c("yearmon","DI")], by="yearmon")
    }else{
      full=data.frame(matrix(ncol=9, nrow=0))
      res=data.frame(matrix(ncol=12, nrow=0))
    }
    
    ### Dispersal Kernel outputs
    if(!is.na(data$release_latitude[1])){
      message('- Calculating dispersal distances and bearings')
      pts<-data.frame(lat=data$latitude, lon=data$longitude); coordinates(pts)<-~lon+lat; projection(pts)<-ll
      pt<-data.frame(lat=data$release_latitude[1], lon=data$release_longitude[1]); coordinates(pt)<-~lon+lat; projection(pt)<-ll
      disp<-data.frame(tag_id=data$tag_id[1], transmitter_id=data$transmitter_id[1], species=data$scientific_name, installation_name=data$installation_name, station_name=data$station_name, 
                          ReleaseDate=as.POSIXct(strptime(data$ReleaseDate, format="%Y-%m-%d %H:%M:%S"), tz="GMT"),
                          ReleaseLat=data$release_latitude, ReleaseLon=data$release_longitude,
                          detection_timestamp=data$dt, lat=data$latitude, lon=data$longitude)
      ### Straight line distance between release location and each detection
      disp$disrel<-spDistsN1(pts=pts, pt=pt,longlat=TRUE)
      disp$azrel<-sapply(1:nrow(disp), function(x) azrel=gzAzimuth(from=matrix(c(disp$ReleaseLon[x],disp$ReleaseLat[x]),1,2), to=matrix(c(disp$lon[x],disp$lat[x]),1,2)))
      disp$azrel<- ifelse(!is.na(disp$azrel) & disp$azrel<0, disp$azrel+360, disp$azrel)
      ### Straight line distance between consecutive detections
      disp$discon<-c(spDistsN1(matrix(c(disp$ReleaseLon[1],disp$ReleaseLat[1]),1,2), matrix(c(disp$lon[1],disp$lat[1]),1,2), longlat=TRUE),
                        sapply(2:nrow(disp), function(x) discon=spDistsN1(matrix(c(disp$lon[x-1],disp$lat[x-1]),1,2), matrix(c(disp$lon[x],disp$lat[x]),1,2), longlat=TRUE)))
      disp$azcon<-c(gzAzimuth(from=matrix(c(disp$ReleaseLon[1],disp$ReleaseLat[1]),1,2), to=matrix(c(disp$lon[1],disp$lat[1]),1,2)),
                       sapply(2:nrow(disp), function(x) azcon=gzAzimuth(from=matrix(c(disp$lon[x-1],disp$lat[x-1]),1,2), to=matrix(c(disp$lon[x],disp$lat[x]),1,2))))
      disp$azcon<- ifelse(!is.na(disp$azcon) & disp$azcon<0, disp$azcon+360, disp$azcon)
    }else{
      message('- Calculating dispersal distances and bearings')
      pts<-data.frame(lat=data$latitude, lon=data$longitude); coordinates(pts)<-~lon+lat; projection(pts)<-ll
      disp<-data.frame(tag_id=data$tag_id[1], transmitter_id=data$transmitter_id[1], species=data$scientific_name, installation_name=data$installation_name, station_name=data$station_name, 
                          ReleaseDate=NA, ReleaseLat=NA, ReleaseLon=NA,
                          detection_timestamp=data$dt, lat=data$latitude, lon=data$longitude)
      ### Straight line distance between release location and each detection
      disp$disrel<-NA; disp$azrel<-NA
      ### Straight line distance between consecutive detections
      disp$discon<-c(NA, sapply(2:nrow(disp), function(x) discon=spDistsN1(matrix(c(disp$lon[x-1],disp$lat[x-1]),1,2), matrix(c(disp$lon[x],disp$lat[x]),1,2), longlat=TRUE)))
      disp$azcon<-c(NA, sapply(2:nrow(disp), function(x) azcon=gzAzimuth(from=matrix(c(disp$lon[x-1],disp$lat[x-1]),1,2), to=matrix(c(disp$lon[x],disp$lat[x]),1,2))))
      disp$azcon<- ifelse(!is.na(disp$azcon) & disp$azcon<0, disp$azcon+360, disp$azcon)
    }
    
    ## Calculate and add activity space metrics (MCP and BBKUD) for each yearmon
    ## Calculate centre of activity (COA)
    message('- Calculating Center of Activity positions')
    step<-timestep*60 ## converts timestep from minutes to seconds
    data$DateTime<-cut(data$dt, breaks=seq(from=trunc(min(data$dt, na.rm=TRUE), "day"), to=trunc(max(data$dt, na.rm=TRUE), "day")+86400, by=step))
    cenac<-ddply(data, "DateTime", summarize, 
                 tag_id=NA, transmitter_id=NA,
                 yearmon=yearmon[1],
                 species=scientific_name[1],
                 meanlat=mean(latitude), meanlon=mean(longitude))
    cenac$tag_id<-data$tag_id[1]; cenac$transmitter_id<-data$transmitter_id[1]
    cenac<-cenac[!is.na(cenac$meanlat),]
    
    ### remove subsets with fewer than 5 detections
    COA<-droplevels(cenac[cenac$yearmon%in%names(table(cenac$yearmon)[table(cenac$yearmon)>5]),])
    cenac$yearmon<-NULL
    
    if(nrow(unique(COA[,c("meanlat","meanlon")]))>5){
      message('- Estimating overall and subsetted MCP and BBKUD')
      ## Setup spatial data and convert from lat long to UTM
      sdat<-COA
      coordinates(sdat)<-~meanlon+meanlat; projection(sdat)<-ll
      dat<-spTransform(sdat,utm)
      
      ### subset data for each yearmon and calculate MCParea and bbkud areas
      ### MCP area calculation. Adds MCP area values into res data frame
      mcparea<-mcp(dat[,"yearmon"], percent=100, unin="m", unout="m2")@data; colnames(mcparea)<-c("yearmon","mcp")
      res<-merge(res,mcparea, by="yearmon", all=TRUE)
      mcputm<-mcp(dat[,"tag_id"], percent=100, unin="m", unout="m2")
      mcpcont<-spTransform(mcputm, ll)
      full$mcp<-mcputm@data$area
      
      ###### Brownian Bridge KUD area calculation
      ### Define grid
      width<-ceiling(max((extent(dat)[2]-extent(dat)[1])/2,(extent(dat)[4]-extent(dat)[3])/2))
      xcen<-(extent(dat)[2]+extent(dat)[1])/2; ycen<-(extent(dat)[4]+extent(dat)[3])/2
      gr<-expand.grid(x=seq(xcen-(width*ext), xcen+(width*ext),len=grid),y=seq(ycen-(width*ext), ycen+(width*ext),len=grid)); coordinates(gr) <- ~x+y; gridded(gr)<-TRUE 
      
      ### Full data
      tf<-as.ltraj(xy=coordinates(dat), date=ymd_hms(dat$DateTime, tz="GMT"), id=dat$transmitter_id, typeII=TRUE, proj4string=utm)
      s1f<-(liker(tf, rangesig1=c(0,500), sig2=sig2, byburst=FALSE, plotit=FALSE)[[1]]$sig1)/div
      tryCatch({
        kbfull<-kernelbb(tf, sig1=s1f, sig2=sig2, grid=gr)
        bf<-kernel.area(kbfull, percent=c(50,95), unin="m", unout="m2")
        full$bbk50<-bf[1]; full$bbk95<-bf[2] 
        },error=function(e){message("ERROR in calculating full BBKUD estimates and area:",conditionMessage(e))})
      
      ### Subsetted data
      ym<-ddply(COA, "yearmon", function(x) nrow(unique(x[,c("meanlat","meanlon")]))); yfac<-ym[ym$V1>5,"yearmon"]
      if(length(yfac)>0){
        traj<-as.ltraj(xy=coordinates(dat), date=ymd_hms(dat$DateTime, tz="GMT"), id=dat$yearmon, typeII=TRUE, proj4string=utm)[yfac]
        s1<-liker(traj, rangesig1=c(0,500), sig2=sig2, byburst=FALSE, plotit=FALSE)
        sig1<-(unname(sapply(s1, '[[', 1)))/div
        bb<-NA
        tryCatch({
          kbb<-kernelbb(traj, sig1=sig1, sig2=sig2, grid=gr)
          bb<-kernel.area(kbb, percent=c(50,95), unin="m", unout="m2")
          },error=function(e){message("ERROR in calculating subsetted BBKUD estimates and area:",conditionMessage(e))})
        barea<-data.frame(yearmon=as.character(ym[ym$V1>5,"yearmon"]), bbk50=unlist(data.frame(bb)[1,]), bbk95=unlist(data.frame(bb)[2,]))
        res<-merge(res,barea, by="yearmon", all=TRUE)
      }else{
        message('WARNING: Not enough unique relocations to estimate subsetted BBKUDs')
        res$bbk50<-NA; res$bbk95<-NA
      }
      
      ##### Cumulative MCP and KUDs
      ##### Cumulative trajectaries for cumulative MCP and KUDs
      if(cumulative){
        message('- Estimating cumulative activity space metrics:')
        cumareas<-data.frame(yearmon=levels(dat$yearmon), cmcp=NA, ck50=NA, ck95=NA)
        for(c in 1:length(levels(dat$yearmon))){
          cdat<-subset(dat, yearmon%in%levels(dat$yearmon)[1:c])
          cumareas[c,"cmcp"]<-mcp.area(cdat[,"transmitter_id"],percent=100, unin="m", unout="m2", plotit=FALSE)
          ct<-as.ltraj(xy=cdat@coords, date=ymd_hms(cdat$DateTime, tz="GMT"), id=cdat$transmitter_id, typeII=TRUE, proj4string=utm)
          csig1<-(liker(ct, rangesig1=c(0,500), sig2=sig2, byburst=FALSE, plotit=FALSE)[[1]]$sig1)/div
          tryCatch({
            ckbb<-kernelbb(ct, sig1=csig1, sig2=sig2, grid=gr)
            cumareas[c,c("ck50","ck95")]<-kernel.area(ckbb, percent=c(50,95), unin="m", unout="m2")
            },error=function(e){message("ERROR in calculating cumulative BBKUD estimates and area:",conditionMessage(e))})
          setTxtProgressBar(txtProgressBar(min=0, max=length(levels(dat$yearmon)), style=3), c)
        }
        res<-merge(res, cumareas, by="yearmon", all=TRUE)
      }
      
      if(plotfull){
        message('- Plotting full BBKUD estimates with COA points')
        dev.new(noRStudioGD=TRUE)
        plot(gr, col=NA)
        plot(raster(getvolumeUD(kbfull)), col=colorRampPalette(brewer.pal('Spectral', n=11))(400), zlim=c(0,100))
        tryCatch({plot(getverticeshr(kbfull, 95), border=1, lty=3, add=TRUE)},error=function(e){message("ERROR in plotting 95% overall BBKUD contour:",conditionMessage(e))})
        tryCatch({plot(getverticeshr(kbfull, 50), border=1, lty=1, add=TRUE)},error=function(e){message("ERROR in plotting 50% overall BBKUD contour:",conditionMessage(e))})
        points(dat, pch=20, cex=0.5, col="blue")
        plot(mcputm, col=NA, border=3, add=TRUE)
        map.axes()
        legend('bottomleft', pch=c(NA,NA,NA,20), col=c(1,1,'green','blue'), lty=c(3,1,1,NA), 
               legend=c("BBKUD95%","BBKUD50%","MCP","COA"), cex=0.8, bg="white")
      }
      
      if(plotsub){
        message('- Plotting subsetted KUD estimates')
        if(length(yfac)>1){
          ras<-stack(lapply(getvolumeUD(kbb), raster))
          pdat<-lapply(1:length(traj), function(x) traj[[x]][,c(1,2)])
          tryCatch({sub95<-getverticeshr(kbb, 95)},error=function(e){message("ERROR in plotting 95% subsetted bbKUD contour:",conditionMessage(e))})
          tryCatch({sub50<-getverticeshr(kbb, 50)},error=function(e){message("ERROR in plotting 50% subsetted bbKUD contour:",conditionMessage(e))})
          dev.new(noRStudioGD=TRUE)
          plot(rasterVis::levelplot(ras,par.settings=rasterTheme(region=colorRampPalette(brewer.pal('Spectral', n=11))(100)), contour=FALSE, colorkey=TRUE, names.attr=as.character(yfac),
                                    panel=function(...){ 
                                      panel.levelplot(...) 
                                      if(exists("sub95")){for(s in 1:length(sub95@polygons[[panel.number()]]@Polygons)){panel.polygon(sub95@polygons[[panel.number()]]@Polygons[[s]]@coords, border=1, lty=3)}}
                                      if(exists("sub50")){for(s in 1:length(sub50@polygons[[panel.number()]]@Polygons)){panel.polygon(sub50@polygons[[panel.number()]]@Polygons[[s]]@coords, border=1, lty=1)}}
                                      panel.points(pdat[[panel.number()]], pch=20, cex=0.1, col="blue")}))
        }else{
            message("   ERROR: not enough temporal subsets for a subsetted plot")
        }
      }
      
      if(storepoly){
        message('- Shapefiles created and stored in output')
        out<-list(full=full,subset=res, COA=cenac, disp=disp, sp=list(raster_full=NA,raster_sub=NA,mcpcont=mcpcont))
        tryCatch({ras_full<-raster(getvolumeUD(kbfull)); projection(ras_full)<-utm; out$sp$raster_full<-projectRaster(ras_full, crs=ll) },error=function(e){message("ERROR in saving full BBKUD raster:",conditionMessage(e))})
        tryCatch({ras_sub<-ras<-stack(lapply(getvolumeUD(kbb), raster)); projection(ras_sub)<-utm; out$sp$raster_sub<-projectRaster(ras_sub, crs=ll) },error=function(e){message("ERROR in saving subsetted BBKUD rasterstack:",conditionMessage(e))})
        }else{
          out<-list(full=full, subset=res, COA=cenac, disp=disp)
          }
      }else{
        message('WARNING: Less than 5 unique relocations, cannot estimate MCP and BBKUDs')
        full$mcp<-NA; full$bbk50<-NA; full$bbk95<-NA
        res$mcp<-NA; res$bbk50<-NA; res$bbk95<-NA
        out<-list(full=full, subset=res, COA=cenac, disp=disp)
        }
    }else{
      message('ERROR: Not enough reliable relocation data to calculte movement and activity space metrics\n- Check tag details are included in taginfo file\n- Check Detection QC values of IMOS data')
      out<-list(full=data.frame(matrix(ncol=12, nrow=0)),
                subset=data.frame(matrix(ncol=14, nrow=0)),
                COA=data.frame(matrix(ncol=6, nrow=0)),
                disp=data.frame(matrix(ncol=15, nrow=0)))
      }
  return(out)
}

