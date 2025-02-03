# app.R
library(shiny)
library(ggplot2)
library(data.table)

# UI definition

sam_data <- readRDS("sam-cmip6.Rds")

source("plot_sam.R")

ui <- fluidPage(
  titlePanel("Trend Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      # Multiple selection for seasons
      checkboxGroupInput("seasons",
                         "Select Seasons:",
                         choices = c("DJF", "MAM", "JJA", "SON"),
                         selected = c("DJF", "MAM", "JJA", "SON")),
      
      # Slider for cut year
      sliderInput("cut_year",
                  "Select Cut Year:",
                  min = 1979,
                  max = 2014,
                  value = 1999,
                  step = 1),
      
      # Checkbox for common models
      checkboxInput("common_models",
                    "Only models with all forcings?",
                    value = TRUE),
      
      # Add a plot button to prevent excessive recomputation
      actionButton("update_plot", "Update Plot")
    ),
    
    mainPanel(
      # Plot output
      plotOutput("trendPlot")
    )
  )
)

# Server logic
server <- function(input, output) {
  
  # Reactive expression for the plot
  plot_data <- eventReactive(input$update_plot, {
    # Validate that at least one season is selected
    validate(
      need(length(input$seasons) > 0, "Please select at least one season")
    )
    
    # Call your plotting function with the user inputs
    plot_sam(data = sam_data,  # You'll need to load your data
             seasons = input$seasons,
             cut_year = input$cut_year,
             common_models = input$common_models)
  })
  
  # Render the plot
  output$trendPlot <- renderPlot({
    plot_data()
  })
}

# Run the app
shinyApp(ui = ui, server = server)