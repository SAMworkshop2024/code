# app.R
library(shiny)
library(ggplot2)
library(data.table)

# UI definition


source("plot_sam.R")

ui <- fluidPage(
  titlePanel("Trend Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      # Multiple selection for seasons
      checkboxGroupInput("seasons",
                         "Select Seasons:",
                         choices = c("DJF", "MAM", "JJA", "SON"),
                         selected = c("DJF", "JJA")),
      
      checkboxGroupInput("forcing",
                         "Select forcings:",
                         choices = c("historical", "hist-GHG", "hist-stratO3", "hist-totalO3"),
                         selected = c("historical", "hist-GHG", "hist-stratO3")),
      
      # Slider for cut year
      sliderInput("cut_year",
                  "Select Cut Year:",
                  min = 1979,
                  max = 2014,
                  value = 1999,
                  step = 1),
      
      # Checkbox for common models
      checkboxInput("common_models",
                    "Only models with all forcings",
                    value = TRUE),
      checkboxInput("continuous",
                    "Continuos regression",
                    value = TRUE)
    ),
    
    mainPanel(
      # Plot output
      plotOutput("trendPlot"),
      plotOutput("linePlot")
    )
  )
)

# Server logic
server <- function(input, output) {
  sam_data <- readRDS("sam-cmip6.Rds")[year(time) %between% c(1979, 2014)] |> 
    _[, .(sam = mean(sam)), by = .(time = seasonally(time),
                                   institute, model, member, forcing)]
  sam_data[, season := season(time)]
  
  era5 <- sam_data[model == "ERA5"]
  sam_data <- sam_data[model != "ERA5"]
  
  sam_data[, forcing := factor(forcing, 
                               levels = c("historical", "hist-GHG", "hist-stratO3", "hist-totalO3"))]
  
  data_experiments <- reactive(sam_data[forcing %in% input$forcing])
  
  data_models <- reactive(select_models(data_experiments(), input$common_models))
  data_season <- reactive(select_season(data_models(), input$seasons))
  era5_season <- reactive(select_season(era5, input$seasons))
  
  # Reactive expression for the plot
  plot_data <- reactive({
    plot_sam(data = data_season(),  # You'll need to load your data
             obs = era5_season(),
             cut_year = input$cut_year,
             continuous = input$continuous)
  })

  # Render the plot
  output$trendPlot <- renderPlot({
    plot_data()[[1]]
  })
  
  output$linePlot <- renderPlot({
    plot_data()[[2]]
  })
  
  
}

# Run the app
shinyApp(ui = ui, server = server)

