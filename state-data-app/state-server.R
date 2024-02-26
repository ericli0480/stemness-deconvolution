
server <- function(input, output, session) {
  
  # Allow user to choose column of dataset
  output$datasetcolumn2 <- renderUI({
    chosen_dataset <- read.csv(paste("data-oh/", input$datasetname, sep = ""))
    selectInput("datasetcolumn", "Data choice:",
                # Only allow user to plot numerical data
                Filter(function(x) class(chosen_dataset[,match(x, names(chosen_dataset))])=="numeric",
                       names(chosen_dataset))
    )
  })
  
  
  
  # Map output for Ohio
  output$map_oh <- renderLeaflet({
    
    # Read shapefiles
    load("shapefile_oh_pa.RData")

    chosen_dataset_oh <- read.csv(paste("data-oh/", input$datasetname, sep = ""))
    chosen_column <- input$datasetcolumn
    chosen_index <- match(chosen_column, names(chosen_dataset_oh))

    # Merge shapefile with desired data
    leaflet_oh <- merge(x = shapefile_oh,
                        y = chosen_dataset_oh,
                        by.x = "geoid", by.y = "FIPS",
                        all.x = TRUE)

    # Color palette for shading
    try(
      pal <- colorNumeric(palette = "Blues", domain = chosen_dataset_oh[,chosen_index]),
      silent = TRUE # Suppress warning from data not loading yet
    )
    
    try ({ # Suppress warning from data not loading yet
      names(leaflet_oh)[which(names(leaflet_oh) == chosen_column)] <- "X"
      leaflet() %>%
        addTiles() %>%
        addPolygons(data = leaflet_oh,
                    stroke = FALSE,
                    smoothFactor = 0.2,
                    fillOpacity = 1,
                    color = ~pal(X),
                    label = ~paste0(namelsad, "=", X)
        )
      }, silent = TRUE)
  })
  
  
  
  # Map output for Pennsylvania
  output$map_pa <- renderLeaflet({
    
    # Read shapefiles
    load("shapefile_oh_pa.RData")
    
    chosen_dataset_pa <- read.csv(paste("data-pa/", input$datasetname, sep = ""))
    chosen_column <- input$datasetcolumn
    chosen_index <- match(chosen_column, names(chosen_dataset_pa))
    
    # Merge shapefile with desired data
    leaflet_pa <- merge(x = shapefile_pa,
                        y = chosen_dataset_pa,
                        by.x = "geoid", by.y = "FIPS",
                        all.x = TRUE)
    
    # Color palette for shading
    try(
      pal <- colorNumeric(palette = "Blues", domain = chosen_dataset_pa[,chosen_index]),
      silent = TRUE # Suppress warning from data not loading yet
    )
    
    try ({ # Suppress warning from data not loading yet
      names(leaflet_pa)[which(names(leaflet_pa) == chosen_column)] <- "X"
      leaflet() %>%
        addTiles() %>%
        addPolygons(data = leaflet_pa,
                    stroke = FALSE,
                    smoothFactor = 0.2,
                    fillOpacity = 1,
                    color = ~pal(X),
                    label = ~paste0(namelsad, "=", X)
        )
    }, silent = TRUE)
  })
}