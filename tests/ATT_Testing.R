
## Testing script

data(tagdata)
tag1<-dplyr::filter(tagdata, tag_id %in% 51448725 )
data(taginfo)
data(statinfo)

source("R/setupData.R")
source("R/abacusPlot.R")
source("R/detectionSummary.R")
source("R/dispersalSummary.R")
source("R/COA.R")
source("R/HRprocess.R")
source("R/HRSummary.R")

ATTdata<-setupData(tagdata, taginfo, statinfo)
# ATTdata<-setupData(tag1, taginfo, statinfo)

d<-detectionSummary(ATTdata)
dd<-dispersalSummary(ATTdata)
abacusPlot(ATTdata)
COAdata<-COA(ATTdata, split=F)

hr<- HRSummary(COAdata, projCRS=CRS("+init=epsg:3577"), type="BBKUD", storepoly=T)


# prep<- dat %>%
#   group_by(Tag.ID) %>%
#   do(HRprocess(., ll=CRS("+init=epsg:3577"), utm=projCRS))

