#' Script for caching large queries
#'
#' The top level visualizations for the SHAP analysis and input data load
#' a lot of data. The process is very time intensive, therefor this scripts
#' caches all required data in a compact format
#'
#' Author: julian.kreis@external.merkcgroup.com
#' Date: Mon Aug  2 07:44:11 2021
#' --------------
library(dplyr)

# switch to app folder for caching to correct folder
setwd(here::here("inst/app/"))

view = "merck_cell_var_guan_features_21q2_dev_kreis_0.fst"

feature_annotations = create_feature_categories(view)

# iterate through top level categories
categories = feature_annotations$category %>%
  unique()

# call top level data
top_level =  summarize_input_values(feature_annotation = feature_annotations,
                                    view = view,
                                    active_x = NULL)

# iterate through top categories
sec_level =  purrr::map(categories, summarize_input_values, view=view,
                        feature_annotation=feature_annotations)

# call Expression data
expression =  summarize_input_values(feature_annotation = feature_annotations,
                                     view = view,
                                     active_x = "Expression")


shap_view = "merck_cell_var_synergy_efficacy_21q2_dev_kreis_0.fst"
input_view = "merck_cell_var_guan_features_21q2_dev_kreis_0.fst"

shap_call = function(x, readout) {
  summarize_shap_values(feature_annotation = feature_annotations,
                        shap_view = shap_view,
                        input_view = input_view,
                        active_x = x,
                        readout = readout)
}


# call top level data
top_level_aoc =  shap_call(NULL, "aoc")
top_level_bliss =  shap_call(NULL, "bliss")

# iterate through top categories
sec_level_aoc =  purrr::map(categories, shap_call, "aoc")
sec_level_bliss =  purrr::map(categories, shap_call, "bliss")

# call Expression data
expression_aoc =  shap_call("Expression", "aoc")
expression_bliss =  shap_call("Expression", "bliss")

# call Expression data
level_3 = c("Treatment_1_name", "Treatment_2_name", "Treatment_1_moa",
            "Treatment_2_moa")
level_3_aoc =  purrr::map(level_3, shap_call, "aoc")
level_3_bliss =  purrr::map(level_3, shap_call, "bliss")


# switch back to main folder
setwd(here::here())
