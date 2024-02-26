################################################################################
##  Purpose: Visualize Ohio and Pennsylvania county data.
##  ---------------------------------------------------------------------------
##  Platform, Bugs and Side-Effects:
##    -- R x64 4.3.1 under Mac OS X
##    -- Temporary error message when switching datasets:
##       Error in $: $ operator is invalid for atomic vectors
##  ---------------------------------------------------------------------------
##  Author: Eric Li
##          Elkins High School
##          Missouri City, TX 77459
##          Phone: 281-302-9126
##          Email: ericjli0480@gmail.com
################################################################################

rm(list = ls())

# Comment out setwd for deployment
setwd("/Users/ericli/Dropbox/cwru-project/state-data-app")

library(shiny)
library(leaflet)
library(sf)

source("state-ui.R")
source("state-server.R")

shinyApp(ui = ui, server = server) # Run app


