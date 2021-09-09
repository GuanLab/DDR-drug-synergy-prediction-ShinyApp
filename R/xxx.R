#' @import ggplot2 ggiraph dplyr
#' @importFrom bslib bs_theme
#' @importFrom shiny HTML h2 testServer reactiveVal NS column div span icon actionButton reactive moduleServer req validate need observeEvent updateActionButton fluidRow h4 bindCache tagList fileInput downloadLink p uiOutput renderUI incProgress withProgress wellPanel shinyOptions
#' @importFrom lightgbm readRDS.lgb.Booster saveRDS.lgb.Booster
#' @importFrom stats median quantile setNames var
#' @importFrom utils data read.table write.table
#' @importFrom rlang sym syms :=
#' @importFrom shinyWidgets pickerInput updatePickerInput pickerOptions
#' @importFrom viridis scale_fill_viridis
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom shinyjs enable disable hide show useShinyjs hidden disabled
#' @importFrom DT renderDT DTOutput datatable renderDataTable dataTableProxy selectRows JS
#' @importFrom glue glue
#' @importFrom digest digest
#' @importFrom forcats fct_inorder
#' @importFrom tidyr nest unnest separate gather unite
#' @importFrom broom tidy
#' @importFrom shinyjqui jqui_resizable
#' @importFrom stringr str_sub str_to_title
#' @importFrom R.utils seqToIntervals
#' @importFrom plyr alply rbind.fill.matrix
#' @importFrom magrittr %>%
utils::globalVariables(c(".identifier_sample_name", ".metadata_cancer_subtype",
                         ".metadata_cancer_type", ".metadata_moa_1", ".metadata_moa_2",
                         ".metadata_treatment_1", ".metadata_treatment_2", ".prediction_aoc",
                         ".prediction_bliss", "AOC", "Bliss", "Total", "across", "aoc_features",
                         "bliss_features", "category", "chemical_structure", "connectivity_aoc",
                         "connectivity_bliss", "drug_moa_combn", "estimate", "feature", "feature_class",
                         "from", "gene", "genes_aoc", "genes_bliss", "geneset_cnt", "measure", "med",
                         "median_value", "metadata_moa_1", "metadata_moa_2", "metadata_treatment_1",
                         "metadata_treatment_2", "middle", "molecular_biomarker", "n", "network_genes_aoc",
                         "network_genes_bliss", "shap", "subcategory", "targets_aoc",
                         "targets_bliss", "tibble", "to", "tool_text", "x", "cols"))