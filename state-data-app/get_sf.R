setwd("/Users/ericli/Dropbox/CWRU Project/state-data-app")
library(shiny)
library(leaflet)
library(sf)
shapefile <- st_read("us-county-boundaries/us-county-boundaries.shp")
shapefile_oh <- shapefile[which(shapefile$state_name == "Ohio"),]
shapefile_pa <- shapefile[which(shapefile$state_name == "Pennsylvania"),]

save(shapefile_oh, shapefile_pa, file = "shapefile_oh_pa.RData")
