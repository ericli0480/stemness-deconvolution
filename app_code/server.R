################################################################################
##  Purpose: Implement Method 4 in real data analysis using Shiny
##    -- Weighted least squares with constraints and quadratic programming 
##    -- Cross-validation to select h 
##  server.R file to accompany app.R
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

source("M4functions.R") # Method 4 and cross-validation functions

server <- function(input, output, session) {
  
  # Automatically send user to "Cross-validation" tab when cross-validation is
  # done
  observeEvent(input$CV, { 
    updateTabsetPanel(session, inputId = "alltabs", selected = "cvtab")
  })
  
  # Automatically send user to "Histogram" tab when algorithm finishes running
  observeEvent(input$go, {
    updateTabsetPanel(session, inputId = "alltabs", selected = "histtab")
  })
  
  # Generate new cross-validation plot each time the user inputs CV.lo, CV.hi,
  # CV.numpts
  genCV <- eventReactive(input$CV, {
    req(input$data)
    req(input$spatial)
    req(input$CV.lo)
    req(input$CV.hi)
    req(input$CV.numpts)
    
    # Inputs
    load(input$data$datapath)
    spatial.data <- read.csv(input$spatial$datapath, header = F)
    
    CV.lo <- input$CV.lo
    CV.hi <- input$CV.hi
    CV.numpts <- input$CV.numpts
    
    # Reformat input dataset to make it compatible with our functions
    colnames(inputDataset)[length(colnames(inputDataset))] <- "Y.orig"
    celltypes <- colnames(inputDataset)[2:(length(colnames(inputDataset)) - 1)]
    for (i in 2:(length(colnames(inputDataset)) - 1)) {
      colnames(inputDataset)[i] <- paste("W", i - 1, sep = "")
    }
    
    # Read and reformat the spatial data (x and y coordinates of spots)
    spatial.data <- spatial.data[, c(1, 5, 6)] # Only keep the needed info
    names(spatial.data) <- c("barcodes", "x.coord", "y.coord")
    
    # Put the information in one object
    dat1 <- merge(inputDataset, spatial.data, by = "barcodes")
    
    # Scale down coordinates to between 0 and 1.
    # This will make it easier to generate kernel weights.
    dat1$x.coord <- dat1$x.coord / max(dat1$x.coord)
    dat1$y.coord <- dat1$y.coord / max(dat1$y.coord)
    
    N <- dim(dat1)[[1]] # Number of cells in dataset
    K <- length(celltypes) # Number of cell types in dataset
    
    # Instantiate progress indicator
    # Needed due to longer running time
    CVprogress <- Progress$new()
    on.exit(CVprogress$close())
    CVprogress$set(message = "Running cross-validation...", value = 0)
    
    # Run cross-validation with the user-inputted numbers
    CVplotOutput <- CV.plot(lo = CV.lo, 
                            hi = CV.hi, 
                            num.pts = CV.numpts, 
                            dat1, 
                            K,
                            CVprogress)
    
    # Send results to download button
    output$downloadcv <- downloadHandler( # Download histograms
      filename = "crossvalidation.pdf",
      content = function(file) {
        pdf(file) # Start a new PDF file
        par(mfrow = c(1, 1)) # One plot per page
        print(CVplotOutput)
        dev.off()
      }
    )
    
    CVplotOutput
    
  })
  
  
  # When the user presses the "Deconvolute" button, run method 4 on the
  # inputted data
  deconvolute <- eventReactive(input$go, {
    req(input$data)
    req(input$spatial)
    req(input$h)
    
    # Inputs
    load(input$data$datapath)
    spatial.data <- read.csv(input$spatial$datapath, header = F)
    h <- input$h
    
    # Reformat input dataset
    colnames(inputDataset)[length(colnames(inputDataset))] <- "Y.orig"
    celltypes <-
      colnames(inputDataset)[2:(length(colnames(inputDataset)) - 1)]
    for (i in 2:(length(colnames(inputDataset)) - 1)) {
      colnames(inputDataset)[i] <- paste("W", i - 1, sep = "")
    }
    
    
    # Read the spatial data (x and y coordinates of spots)
    spatial.data <- spatial.data[, c(1, 5, 6)] # Only keep the info we need
    names(spatial.data) <- c("barcodes", "x.coord", "y.coord")
    
    # Put the weights, stemnesses, and coordinates in one object
    dat1 <- merge(inputDataset, spatial.data, by = "barcodes")
    
    # Scale down coordinates to between 0 and 1.
    # This will make it easier to generate kernel weights.
    dat1$x.coord <- dat1$x.coord / max(dat1$x.coord)
    dat1$y.coord <- dat1$y.coord / max(dat1$y.coord)
    
    N <- dim(dat1)[[1]] # Number of cells in dataset
    K <- length(celltypes) # Number of cell types in dataset
    
    # Lists containing the histograms and spatial plots to be outputted
    returnHist <- vector("list", K) 
    returnSpatial <- vector("list", K)
    
    # List of length 2 containing the lists returnHist and returnSpatial
    returnAllPlots <- list()
    
    # Instantiate progress indicator
    # Mostly a formality, Method 4 runs quickly
    M4progress <- Progress$new()
    on.exit(M4progress$close())
    M4progress$set(message = "Running deconvolution...", value = 0)
    
    # Run Method 4 with the user-selected h (optimal bandwidth)
    alpha <- method4.run(h = h, dat1 = dat1, K = K, M4progress)$alpha
    colnames(alpha) <- celltypes
    
    # Create histograms of the alphas, each bar a different color
    # Plot proportions, not counts
    for (i in 1:K) {
      df <- data.frame(x = alpha[, i])
      myplot <- ggplot(df, aes(x)) +
        geom_histogram(bins = 20, 
                       fill = colorRampPalette(c("#ADD8E6", "#003366"))(20)) +
        aes(y = stat(count)/sum(stat(count))) + # Convert to proportions
        ylim(0, 1) + # Make all plots have the same scale
        xlab("Stemness") + ylab("Proportion of Cells") +
        labs(title = celltypes[i])
      returnHist[[i]] <- myplot
    }
    
    # Plot stemness of each cell spatially, using same scale as histograms
    for (i in 1:K) {
      df.x <- dat1$y.coord
      df.y <- -1 * dat1$x.coord
      df <- data.frame(x = df.x,
                       y = df.y,
                       stemness = alpha[, i])
      myplot <-
        ggplot(df, aes(x, y)) +
        geom_star(aes(colour = stemness, fill = stemness), 
                  starshape = "hexagon") + # Plot spots as hexagons
        scale_color_gradient(low = "#ADD8E6", high = "#003366") + 
        scale_fill_gradient(low = "#ADD8E6", high = "#003366") + 
        labs(title = celltypes[i]) +
        xlab(NULL) + ylab(NULL) +
        theme( # Remove x and y axis numbering
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
        )
      returnSpatial[[i]] <- myplot
    }
    
    returnAllPlots <- list(returnHist, returnSpatial)
    
    # Send results to download buttons
    output$downloadhist <- downloadHandler( # Download histograms
      filename = "histograms.pdf",
      content = function(file) {
        pdf(file) # Start a new PDF file
        par(mfrow = c(1, 1)) # One plot per page
        for(i in 1:length(returnHist)) { 
          print(returnHist[[i]])
        }
        dev.off()
      }
    )
    output$downloadspatial <- downloadHandler( # Download spatial plots
      filename = "spatialplots.pdf",
      content = function(file) {
        pdf(file, hei = 4, wid = 5) # Start a new PDF file
        par(mfrow = c(1, 1)) # One plot per page
        for(i in 1:length(returnSpatial)) { 
          print(returnSpatial[[i]])
        }
        dev.off()
      }
    )
    output$downloadalpha <- downloadHandler( # Download alpha values
      filename = "alpha.csv",
      content = function(file) {
        write.csv(alpha, file) # Output alpha as a .csv file
      }
    )
    
    returnAllPlots
  })
  
  output$CVplot <- renderPlot(genCV())
  
  
  
  drawFinalPlots <- reactive({
    # Brute-force to plot all 15 cell types, both histogram and spatial. 
    # Cannot use a loop due to this known bug:
    # https://stackoverflow.com/questions/51893084/renderplot-issue-when-rendering-a-list-of-plots
    
    K <- length(deconvolute()[[1]]) # Number of plots
    
    if (1 <= K) {
      output$hist1 <- renderPlot(deconvolute()[[1]][[1]])
      output$spatial1 <-
        renderPlot(deconvolute()[[2]][[1]], height = 280, width = 350)
    }
    if (2 <= K) {
      output$hist2 <- renderPlot(deconvolute()[[1]][[2]])
      output$spatial2 <-
        renderPlot(deconvolute()[[2]][[2]], height = 280, width = 350)
    }
    if (3 <= K) {
      output$hist3 <- renderPlot(deconvolute()[[1]][[3]])
      output$spatial3 <-
        renderPlot(deconvolute()[[2]][[3]], height = 280, width = 350)
    }
    if (4 <= K) {
      output$hist4 <- renderPlot(deconvolute()[[1]][[4]])
      output$spatial4 <-
        renderPlot(deconvolute()[[2]][[4]], height = 280, width = 350)
    }
    if (5 <= K) {
      output$hist5 <- renderPlot(deconvolute()[[1]][[5]])
      output$spatial5 <-
        renderPlot(deconvolute()[[2]][[5]], height = 280, width = 350)
    }
    if (6 <= K) {
      output$hist6 <- renderPlot(deconvolute()[[1]][[6]])
      output$spatial6 <-
        renderPlot(deconvolute()[[2]][[6]], height = 280, width = 350)
    }
    if (7 <= K) {
      output$hist7 <- renderPlot(deconvolute()[[1]][[7]])
      output$spatial7 <-
        renderPlot(deconvolute()[[2]][[7]], height = 280, width = 350)
    }
    if (8 <= K) {
      output$hist8 <- renderPlot(deconvolute()[[1]][[8]])
      output$spatial8 <-
        renderPlot(deconvolute()[[2]][[8]], height = 280, width = 350)
    }
    if (9 <= K) {
      output$hist9 <- renderPlot(deconvolute()[[1]][[9]])
      output$spatial9 <-
        renderPlot(deconvolute()[[2]][[9]], height = 280, width = 350)
    }
    if (10 <= K) {
      output$hist10 <- renderPlot(deconvolute()[[1]][[10]])
      output$spatial10 <-
        renderPlot(deconvolute()[[2]][[10]], height = 280, width = 350)
    }
    if (11 <= K) {
      output$hist11 <- renderPlot(deconvolute()[[1]][[11]])
      output$spatial11 <-
        renderPlot(deconvolute()[[2]][[11]], height = 280, width = 350)
    }
    if (12 <= K) {
      output$hist12 <- renderPlot(deconvolute()[[1]][[12]])
      output$spatial12 <-
        renderPlot(deconvolute()[[2]][[12]], height = 280, width = 350)
    }
    if (13 <= K) {
      output$hist13 <- renderPlot(deconvolute()[[1]][[13]])
      output$spatial13 <-
        renderPlot(deconvolute()[[2]][[13]], height = 280, width = 350)
    }
    if (14 <= K) {
      output$hist14 <- renderPlot(deconvolute()[[1]][[14]])
      output$spatial14 <-
        renderPlot(deconvolute()[[2]][[14]], height = 280, width = 350)
    }
    if (15 <= K) {
      output$hist15 <- renderPlot(deconvolute()[[1]][[15]])
      output$spatial15 <-
        renderPlot(deconvolute()[[2]][[15]], height = 280, width = 350)
    }
    
    " "
  })
  
  output$finalstatus <- renderText(drawFinalPlots())
}
