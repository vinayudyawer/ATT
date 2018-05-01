
Animal Tracking Toolbox Quick Guide
===================================

Background
------------

Passive telemetry studies use detection patterns of a tagged animal
within a fixed array to understand movement patterns, habitat use and
activity space. Raw detection data are typically used to calculate
metrics of detection (i.e. number of detections, number of days
detected, number of receivers tag was detected on, index of residence),
dispersal (e.g. distances and bearings between consecutive detections;
step distances and turning angles, distances and bearings between each
detection and release site) and activity space (e.g. Minimum Convex
Polygon MCP area, Kernel Utilisation Distribution area), however the
techniques and parameters used to calculate these metrics are often
customised to each study making cross-study comparisons unreliable. Here
we provide a tool to enable standardisation of the calculation of these
commonly used metrics and provide an analytical tool to facilitate.

<br>

<img src="images/Fig1b.png"/>
<sub>Figure 1. Visual summary of workflow to calculate standardised metrics using the Animal Tracking Toolbox.</sub>

<br>
<br>

The Animal Tracking Toolbox (ATT) is a collection of functions created
in the R statistical environment (R Development Core Team 2018) that
calculates standardised metrics of dispersal and activity space from
passive telemetry to enable direct comparisons between animals tracked
within the same study and between studies or locations. The functions
uses individual detection data files alongside tag metadata and receiver
station information to calculate a range of standardised movement and
activity space metrics. This toolbox can be used to calculate and
visualise standardised metrics of movement and activity space within and
between species tracked at multiple locations.

<br>
<br>

<img src="images/Fig2.png"/>
<sub>Figure 2. Overall activity space metric plots for multiple species tagged at multiple locations (a) Yellowfin Bream (n=1), (b) Yellowtail Kingfish (n=1), (c) Grey Reef Shark (n=1) and (d) Bull Shark (n=1) output using the ATT. Coloured points represent Centres of Activity (60 min time steps) with darker shapes representing core activity space (50% contour of Brownian bridge kernel utilisation distribution; BBKUD) and lighter shapes representing the extent of activity space (95% contour of BBKUD). Black polygons represent overall Minimum Convex Polygons from detection data. Open circles represent locations of VR2W receivers deployed within the IMOS ATF infrastructure and associated research installations.</sub>

<br>
<br>

The ATT was developed to preprocess and calculate standardised metrics
of dispersal and activity space from large-scale detection data housed in
the Integrated Marine Observing System’s Animal Tracking Facility (IMOS
ATF) national data repository. The ATT accepts detection data (referred
to as ‘tagdata’ in the function) exported from the IMOS ATF database
(can be accessed through the [AODN portal](https://portal.aodn.org.au)),
however we are working on including functionality for data export formats 
from the VEMCO data management software VUE and other passive telemetry 
networks. 

<br>
<br>

This manual will outline the required data formats for input 
‘tagdata’ and associated tag metadata (referred to as ‘taginfo’ 
) and receiver station information ('statinfo'). This manual will
also demonstrate how to run the function for a single tag as well as
running the function for a large number of tags.

Installation
------------

``` r
# Currently the development version can be accessed from GitHub:
install.packages("devtools")
devtools::install_github("vinayudyawer/ATT")
```

Functions within the toolbox
------------

The Animal Tracking Toolbox is comprised of five main functions that work in series:

1.  **`setData()`** sets up data and produces a single list 'ATT' object so detection data, tag metadata and station information are all in one place. Initialises data for use with other functions in the toolbox.

2.  **`detectionSummary()`** calculates standard detection metrics using an 'ATT' object. Produces a list with detection metrics calculated over the full tag life and within user-defined temporal subset (i.e. monthly and weekly metrics).

3.  **`dispersalSummary()`** calculates standard dispersal metrics using an 'ATT' object. Produces a tibble dataframe with dispersal distance and bearing measurements between consecutive detections as well as between each detection and release location (if provided in 'taginfo').

4.  **`COA()`** estimates short-term Centers of Activity using an 'ATT' object. Based on technique described in [Simpfendorfer et al. 2002](http://www.nrcresearchpress.com/doi/abs/10.1139/f01-191#.WuggLS_L2XQ). Produces a 'COA' tibble dataframe object with centers of activity estimated within user-defined timesteps.

5.  **`HRSummary()`** calculates standardised activity space metrics using a 'COA' object. Produces a list with activity space metrics calculated over the full tag life and within user-defined temporal subsets (i.e. monthly and weekly). Technique of calculating activity space metrics include minimum convex polygons (*MCP*), fixed kernel utilisation distributions (*fKUD*) or Brownian bridge kernel utilisation distributions (*BBKUD*). Cumulative metrics of activity space is also calculated with `cumulative` argument. Spatial polygons and raster objects for further plotting are also produced with `storepoly` argument.

<br>

In addition to these functions, there are additional functions to help plot detection summaries using an abacus plot (**`abacusPlot()`**). We are working on more functions to help visualise dispersal summaries and activity spaces calculated... Watch this space!!

<br>
<br>

Usage
------------

Setting up data

```{r, include=TRUE, eval=TRUE}
## Input example datasets
data(tagdata) 
data(taginfo)
data(statinfo)

## Setup data for use with the Animal Tracking Toolbox
ATTdata<- setupData(Tag.Detections = tagdata, Tag.Metadata = taginfo, Station.Information = statinfo)

```

Calculating detection metrics
```{r, include=TRUE, eval=TRUE}
## Calculate detecion metrics with monthly subsets chosen
detSum<-detectionSummary(ATTdata, sub = "%Y-%m")

## Accessing metrics of detection for full tag life
detSum$Overall

## Accessing metrics of detection for each temporal subset
detSum$Subsetted

## Create an abacus plot
abacusPlot(ATTdata)

```
<img src="images/Fig3.png"/>

```{r, include=TRUE, eval=TRUE}
## Create a facetted abacus plot for a single individual




```

Calculating dispersal metrics
```{r, include=TRUE, eval=TRUE}
## Calculate dispersal metrics
dispSum<-dispersalSummary(ATTdata)

## Accessing metrics of dispersal
dispSum

```

Calculating activity space metrics
```{r, include=TRUE, eval=TRUE}
## First, estimate Short-term center of activities
COAdata<-COA(ATTdata)

## HRSummary() requires calculation of COAs first
## Estimate 100% MCP areas
mcp_est<-HRSummary(COAdata, projCRS=CRS("+init=epsg:3577"), type="MCP", cont=100)

##*** Warning: the following might take a while to run! ***##
## Estimate 50% and 95% fKUD areas with cumulative metrics calculated
kud_est<-HRSummary(COAdata, projCRS=CRS("+init=epsg:3577"), type="fKUD", cumulative=TRUE)

## Estimate 20%, 50% and 95% BBKUD contour areas and store polygons
kud_est<-HRSummary(COAdata, projCRS=CRS("+init=epsg:3577"), type="BBKUD", cont=c(20,50,95), storepoly=TRUE)

```
***More functions to visualise standardised metrics coming soon!!***

<br>
<br>
<br>

Current version
---------------

1.0.0 (1 May 2018)
