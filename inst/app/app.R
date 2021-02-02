#' 
#' 
#' TODO:
#' - put in quantile normalization ? TMM/RLE?
#' - select moa chemical 
#' - agrregate across different treatement 
#' - a set of samples
#' - cancer type/tissue/moa; summary graph ; box plots
#' - take new cell lines                                  | done
#' - pick cell lines for specific moa combinations
#' - surrogate model (important genes only)               | done
#' 
#' Author: julian.kreis@external.merkcgroup.com
#' Date: Tue Feb  1 17:00:33 2021
#' --------------


library(synergyddrapp)
library(dplyr)

# Define UI for application that draws a histogram
ui = basicPage(
   # App title
   titlePanel("Predict optimal drug combinations"),
   
   # Sidebar layout for input and output definitions
   sidebarLayout(
      # Sidebar panel for inputs
      sidebarPanel(
         shiny::selectInput("model", "Cell Line", choices = character(0),
                            selected = character(0)),
         
         # Input: select a file
         fileInput("file1", "Upload patient's molecular biomarker readouts (tsv file)",
                   multiple = TRUE,
                   accept = c("text/table/tsv",
                              "text/tab-separated-values,text/plain",
                              ".tsv"))
      ),
    mainPanel(
       shiny::fluidRow(
          shiny::plotOutput("combination_efficacy") 
       ),
       shiny::fluidRow(
          DT::dataTableOutput("mytable")
       )
    )  
      
   )
)
   
# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
   con = shiny::callModule(xopmodules::xopdata_profiler_access,
                           "login",
                           quick_access_functions = FALSE)
   
   shiny::observeEvent(con$login(), {
      
      options(xopdata_cache=TRUE)
      lines = xopdata::load(table    = "depmap_ccle_rna_wts_gexp",
                            cols     = "id_cell_line_model",
                            distinct = TRUE, 
                            order = TRUE) %>%
         dplyr::pull()
      options(xopdata_cache=FALSE)
      
      shiny::updateSelectInput(session = session, "model", 
                               choices = lines, selected = character(0))
   })
   
   input_data = shiny::reactiveVal()
   
   shiny::observeEvent(input$model, {
      shiny::req(input$model != "")
      
      genes = c(genes_aoc, genes_bliss) %>%
         unique() %>%
         gsub("_exp", "", .)
      
      # access data
      exp_data = xopdata::load(table    = "depmap_ccle_rna_wts_gexp",
                               cols     = c("id_cell_line_model",
                                            "quant_tpm", "id_gene_symbol"),
                               filter   = list(id_gene_symbol=genes,
                                               id_cell_line_model=input$model),
                               maturity = "dev",
                               version  = 1) %>%
         dplyr::collect()
      
      # reformat
      exp_data = exp_data %>%
         dplyr::mutate(id_gene_symbol = paste(id_gene_symbol, "exp", sep="_")) %>%
         tidyr::spread(key=id_gene_symbol, value="quant_tpm") %>%
         tibble::column_to_rownames("id_cell_line_model")
      
      input_data(exp_data)
   })
   
   shiny::observeEvent(input$file1,{
      
      shiny::req(input$file1)
      p_data   <- read.table(input$file1$datapath, header = TRUE)
      input_data(p_data)
   })
   
   predictions = shiny::reactiveVal()
   
   shiny::observeEvent(input_data(), {
      shiny::req(!is.null(input_data()))
      df_drugs <- predict_optimal_drug_combination(input_data())
      predictions(df_drugs)
   })
   
   output$combination_efficacy <- shiny::renderPlot({
      
      df = predictions()
      
      shiny::req(!is.null(df))
      
      # access the contribution of each feature (list of tibbles)
      five_fold_contribution = df$aoc_contribution
      
      # identify feature names
      feature_name = colnames(five_fold_contribution[[1]])
      
      # get the checmical features
      chemical_features = grep("(MACS|Morgan|RDK|FP2|FP4)", feature_name,
                                value=T)
      
      # for each fold calculate the mean contribution
      shap = purrr::map(five_fold_contribution, function(x) {
         tmp = x[, chemical_features] %>%
            tidyr::nest(morgan = which(grepl("Morgan", colnames(.))),
                        rdk = which(grepl("RDK", colnames(.))),
                        fp2 = which(grepl("FP2", colnames(.))),
                        fp4 = which(grepl("FP4", colnames(.)))) %>%
            dplyr::mutate_all(~ mean(as.matrix(.x[[1]])))
         }) %>%
         dplyr::bind_rows() %>%
         tidyr::gather(key="chemical", value="mean_shap")
      
      # find the mean contribution over all folds
      shap_mean = shap %>%
         dplyr::group_by(chemical) %>%
         dplyr::summarise(mean=mean(mean_shap)) %>%
         dplyr::arrange(mean)
      
      # order plot
      shap = shap %>%
         dplyr::mutate(chemical=factor(chemical, levels=shap_mean$chemical))
      
      # plot data
      ggplot2::ggplot() +
         ggplot2::geom_boxplot(data = shap, ggplot2::aes(x=chemical, y=mean_shap)) +
         ggplot2::geom_point(data = shap_mean, ggplot2::aes(x=chemical, y=mean),
                             color="red", size=4) +
         ggplot2::coord_flip() +
         ggplot2::theme_bw(base_size = 16)
   })
   
   output$mytable <- DT::renderDataTable({
      
      df = predictions()
      
      shiny::req(!is.null(df))
      DT::datatable(filter = "top", data = df$predictions,
                    options = list(scrollX = TRUE))
      })
}

# Run the application 
shinyApp(ui = ui, server = server)



