
ui <- fluidPage(
  title = "OH PA Datasets",
  titlePanel("Visualizing Ohio Datasets"),
  sidebarPanel(
    selectInput("datasetname", "Select dataset:",
                list.files(paste(getwd(), "/data-oh", sep = "")) # List of datasets
                ),
  ),
  mainPanel(
    uiOutput("datasetcolumn2"), # Dropdown to select data to plot from dataset
    leafletOutput("map_oh"),
    leafletOutput("map_pa")
  )
)