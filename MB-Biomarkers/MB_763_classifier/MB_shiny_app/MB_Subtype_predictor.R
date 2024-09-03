# Load necessary libraries
library(shiny)
library(randomForest)
library(dplyr)
library(readr)
library(tibble)
library(shinycssloaders)

# Increase the file upload size limit (e.g., to 1000 MB)
options(shiny.maxRequestSize = 1000 * 1024^2)

# Define UI
ui <- fluidPage(
  titlePanel("MB Subtype Predictor"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose Methylation Data File (CSV format)", accept = ".csv"),
      actionButton("predict", "Predict Subtype"),
      br(),
      p("This Shiny app predicts 12 medulloblastoma subtypes using a Random Forest model trained on a methylation dataset (n=763) from Cancer Cell (2017 Jun 12;31(6):737-754.e6). To use this app for medulloblastoma subtype prediction, please upload your samples' methylation beta values in CSV format, with features in rows and sample names in columns. The file size should be less than 50 MB. For guidance on formatting your data file, please refer to the example dataset provided."),
      downloadButton("downloadData", "Download Example Dataset")
    ),
    
    mainPanel(
      withSpinner(tableOutput("predictions"))
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Load pre-trained Random Forest model
  rf_model <- readRDS("rf_model.rds")
  
  # Define reactive expression to load and process the uploaded data
  uploaded_data <- reactive({
    req(input$file1)
    methyl_data <- read.csv(input$file1$datapath, row.names = 1)
    
    # Ensure methyl_data is a data frame
    methyl_data <- as.data.frame(methyl_data)
    
    # Match the row names of methyl_data with the row names of rf_model$importance
    model_features <- rownames(rf_model$importance)
    common_features <- intersect(model_features, rownames(methyl_data))
    methyl_data <- methyl_data[common_features, , drop = FALSE]
    
    # Transpose the data so that features become columns and samples become rows
    methyl_data <- as.data.frame(t(methyl_data))
    
    # Add missing columns with zero values
    missing_features <- setdiff(model_features, colnames(methyl_data))
    for (feature in missing_features) {
      methyl_data[[feature]] <- 0
    }
    
    # Reorder columns to match model features
    methyl_data <- methyl_data[, model_features, drop = FALSE]
    
    methyl_data
  })
  
  # Define reactive expression to make predictions
  predictions <- eventReactive(input$predict, {
    req(uploaded_data())
    predict(rf_model, uploaded_data())
  })
  
  # Display the predictions
  output$predictions <- renderTable({
    req(predictions())
    data.frame(Sample_ID = rownames(uploaded_data()), Predicted_Subtype = predictions())
  })
  
  # Provide example dataset for download
  output$downloadData <- downloadHandler(
    filename = function() {
      "example_dataset.csv"
    },
    content = function(file) {
      example_data <- data.frame(
        Sample1 = rnorm(10),
        Sample2 = rnorm(10)
      )
      rownames(example_data) <- paste0("Feature", 1:10)
      write.csv(example_data, file)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)
