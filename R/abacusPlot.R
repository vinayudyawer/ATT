
abacusPlot<-function(ATTdata, theme="theme_linedraw", xlab="Date", ylab="Tag ID", det.col=2, tag.col=8, facet=FALSE, new.window=TRUE, ...){

  suppressPackageStartupMessages({require(dplyr); require(ggplot2)})

  if(!inherits(ATTdata, "ATT"))
    stop("Oops! Input data needs to be an 'ATT' object.
         \nSet up your data first using setupData() before running this operation")

  ## Combine Tag.Detection and Tag.Metadata into a combined tibble for plotting
  combdata<- left_join(ATTdata$Tag.Detections, ATTdata$Tag.Metadata, by="Tag.ID")

  ## Find start and end date of taglife
  ss<-combdata %>%
    group_by(Tag.ID) %>%
    summarize(Start = min(first(Release.Date), min(date(Date.Time)), na.rm=T),
              End = max((first(Release.Date)+first(Tag.Life)), max(date(Date.Time)), na.rm=T))

  if(new.window){dev.new(noRStudioGD=TRUE, width=9, height=6)}

  if(facet){
    ggplot(combdata) +
      xlab(xlab) + ylab(ylab) +
      geom_point(aes(x = date(Date.Time), y = as.factor(Station.Name)), col=det.col, ...) +
      facet_wrap(~Tag.ID) +
      scale_x_date(date_labels= "%b\n%Y", minor_breaks = NULL) +
      eval(call(theme))
  }else{
    ggplot(combdata) +
      xlab(xlab) + ylab(ylab) +
      geom_point(aes(x = date(Date.Time), y = as.factor(Tag.ID)), col=det.col, ...) +
      geom_point(data= ss, aes(x = Start, y = as.factor(Tag.ID)), pch="|", col=tag.col, cex=3) +
      geom_point(data= ss, aes(x = End, y = as.factor(Tag.ID)), pch="|", col=tag.col, cex=3) +
      scale_x_date(date_labels= "%b\n%Y", minor_breaks = NULL) +
      eval(call(theme))
  }
}