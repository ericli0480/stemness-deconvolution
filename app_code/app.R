################################################################################
##  Purpose: Implement Method 4 in real data analysis using Shiny
##    -- Weighted least squares with constraints and quadratic programming 
##    -- Cross-validation to select h
##
##  ---------------------------------------------------------------------------
##  Platform, Bugs and Side-Effects:
##    -- Version 18
##    -- R x64 4.3.1 under Mac OS X
##    -- May get the error that matrix D is not positive definite. This is an 
##       issue with the numerical values in the dataset. In these cases our
##       model will not work without adjusting these values slightly.
##  ---------------------------------------------------------------------------
##  Author: Eric Li
##          Elkins High School
##          Missouri City, TX 77459
##          Phone: 281-302-9126
##          Email: ericjli0480@gmail.com
################################################################################

rm(list = ls())

# Shiny app libraries
library(shiny)
library(shinyFiles)
library(bslib)
library(thematic)

# Method 4 libraries
library(ggplot2)
library(quadprog)
library(ggstar)

source("M4functions.R") # Method 4 and cross-validation functions
source("ui.R")
source("server.R")


set.seed(500)
options(shiny.maxRequestSize = 500 * 1024 ^ 2) # Allow larger file uploads



shinyApp(ui = ui, server = server) # Run app