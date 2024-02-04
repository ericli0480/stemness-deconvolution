################################################################################
##  Purpose: Implement Method 4 in real data analysis using Shiny
##    -- Weighted least squares with constraints and quadratic programming 
##    -- Cross-validation to select h 
##  ui.R file to accompany app.R
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

# Shiny app libraries
library(shiny)
library(shinyFiles)
library(bslib)
library(thematic)

# Method 4 libraries
library(ggplot2)
library(quadprog)
library(ggstar)

ui <- fluidPage(
  title = "MyTitle",
  theme = bs_theme(version = 5, bootswatch = "cerulean"),
  
  
  titlePanel(
    "A Deconvolution Model to Estimate Cell-type Specific Stemness from Spatial Transcriptomics Data"
  ),
  
  
  sidebarLayout(
    position = "left",
    
    # Sidebar for user to upload data
    sidebarPanel(
      fileInput("data", "Upload Dataset [.RData]", accept = ".RData", 
                buttonLabel = icon("folder-open")),
      fileInput("spatial", "Upload Spatial Data [.csv]", accept = ".csv",
                buttonLabel = icon("folder-open")),
      numericInput("CV.lo", 
                   label = "Cross-validation: Lowest h", value = NULL),
      numericInput("CV.hi", 
                   label = "Cross-validation: Highest h", value = NULL),
      numericInput("CV.numpts", 
                   label = "Cross-validation: Number of Points", value = NULL),
      actionButton("CV", 
                   label = list(icon("play"), "Cross-validate")),
      numericInput("h", 
                   label = "Selected value of h", value = NULL),
      actionButton("go", 
                   label = list(icon("play"), "Deconvolute")),
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs", id = "alltabs",
                  
                  # "About" tab
                  # Contains brief description of the program and specifies input format
                  tabPanel("About", id = "abouttab",
                           br(),
                           strong("About this application"),
                           p(
                             "The metric of \"stemness\" measures the level of development of cells. In cancer research, it is often valuable to analyze the stemnesses of cells in tissue in order to better understand how various cell types behave in terms of differentiation. Stemness can be calculated using the CytoTRACE algorithm, developed by Stanford University."
                           ),
                           p(
                             "While CytoTRACE was designed to be used on single-cell experiments, it can also be used on spatial transcriptomics datasets, which supplement single-cell experiments with spatial data at the cost of resolution. Thus when using CytoTRACE to spatially analyze stemness, the stemnesses we find are not the stemnesses of every cells, but the average stemnesses of groups of cells, called spots."
                           ),
                           p(
                             "In this application, we apply a deconvolution model to estimate the individual stemnesses of each cell in the dataset, under the reasonable assumption that cells of the same type at the same spot have the same stemness."
                           ),
                           br(),
                           strong("Input format: "),
                           p(
                             "The input dataset should be a .RData file containing a single data frame inputDataset."
                           ),
                           p(
                             "    - The first column should be called \"barcodes\" and contain the DNA barcodes of each spot in the dataset."
                           ),
                           p(
                             "    - The next K columns should contain the cell type proportions for each spot, where K is the number of distinct cell types. These proportions can be found from the CARD algorithm. Each column should be labeled with the name of the cell type. Up to K = 15 cell types are supported."
                           ),
                           p(
                             "    - The final column should be called \"stemness\" and contain the stemness of each spot, found from the CytoTRACE algorithm."
                           ),
                           p(
                             "    - The spatial data is the tissue_positions_list.csv file found in the Visium dataset folder."
                           )
                  ),
                  
                  # "Cross-validation" tab
                  # Contains the outputted cross-validation plot (U-shaped curve)
                  tabPanel("Cross-validation", value = "cvtab",
                           br(),
                           h4("Cross-validation Plot"),
                           plotOutput("CVplot")
                  ),
                  
                  # "Histogram" tab
                  # Contains histogram plots of every cell type's stemness
                  tabPanel("Histograms", value = "histtab",
                           br(),
                           h4("Histograms of Stemness of Each Cell Type"),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%", # 3 plots per row
                               plotOutput("hist1"), plotOutput("hist2"), plotOutput("hist3")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("hist4"), plotOutput("hist5"), plotOutput("hist6")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("hist7"), plotOutput("hist8"), plotOutput("hist9")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("hist10"), plotOutput("hist11"), plotOutput("hist12")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("hist13"), plotOutput("hist14"), plotOutput("hist15")
                             )
                           )
                  ),
                  
                  # "Spatial" tab
                  # Contains spatial plots of stemnesses of every cell type
                  tabPanel("Spatial", value = "spatialtab",
                           br(),
                           h4("Spatial Distribution of Stemness by Cell Type"),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("spatial1"), plotOutput("spatial2"), plotOutput("spatial3")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("spatial4"), plotOutput("spatial5"), plotOutput("spatial6")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("spatial7"), plotOutput("spatial8"), plotOutput("spatial9")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("spatial10"), plotOutput("spatial11"), plotOutput("spatial12")
                             )
                           ),
                           fluidRow(
                             splitLayout(
                               cellWidths = "33%",
                               plotOutput("spatial13"), plotOutput("spatial14"), plotOutput("spatial15")
                             )
                           )
                  ),
                  
                  # "Download" tab
                  # Contains two buttons to download CV plot and histogram +
                  # spatial plots, as well as numerical results alpha
                  tabPanel("Download", value = "downloadtab",
                           br(),
                           h4("Download Results"),
                           downloadButton("downloadcv", 
                                          label = "Cross-validation Plot"),
                           br(), br(),
                           downloadButton("downloadhist", 
                                          label = "Histograms"),
                           br(), br(),
                           downloadButton("downloadspatial", 
                                          label = "Spatial Plots"),
                           br(), br(),
                           downloadButton("downloadalpha", 
                                          label = "Numerical results")
                  )
      )
    )
  ),
  
  verbatimTextOutput("finalstatus") # Indicates that running has finished
)