# example input: 
# molecular biomarkers of treated cell lines



#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
###### call lightGBM in R 
library(lightgbm)

###### call python in R package
library(reticulate)
library(DT)

# load pre-trained lightGBM model


#py_install("pandas") #install pandas in R
#py_install("lightgbm")
#py_install("shap")
reticulate::source_python(here::here("./dependency/main_prediction.py"))
p_data <- read.table(here::here('./dependency/input/example_input.tsv'), header = TRUE)
#df_drugs <- predict_optimal_drug_combination(p_data)

# Define UI for application that draws a histogram
ui = basicPage(
   # App title
   titlePanel("Predict optimal drug combinations"),
   
   # Sidebar layout for input and output definitions
   sidebarLayout(
      # Sidebar panel for inputs
      sidebarPanel(
         # Input: select a file
         fileInput("file1", "Upload patient's molecular biomarker readouts (tsv file)",
                   multiple = TRUE,
                   accept = c("text/table/tsv",
                              "text/tab-separated-values,text/plain",
                              ".tsv"))
      ),
    mainPanel(
       DT::dataTableOutput("mytable")
    )  
      
   )
)
   
# Define server logic required to draw a histogram
server <- function(input, output) {
   output$mytable <- DT::renderDataTable({
      
      shiny::req(input$file1)
      p_data <- read.table(input$file1$datapath, header = TRUE)
      df_drugs <- predict_optimal_drug_combination(p_data)
      return(df_drugs)
      }
                                   )
}

# Run the application 
shinyApp(ui = ui, server = server)

# put in quantile normalization 
# select moa chemical 
 #

# grregate across different treatement 


# a set of samples
# cancer type/tissue/moa; summary graph ; box plots


# take new cell lines

# pick cell lines for specific moa combinations

#surrogate model (important genes only)


