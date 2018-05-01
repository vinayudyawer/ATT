
Animal Tracking Toolbox Quick Guide
===================================

Background
--------

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
<br>

 <div style="text-align:center"><img src="images/Fig1b.png" width="700" />

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

<img src="images/Fig2.png" width="600" /></div>

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

Input data format
------------




