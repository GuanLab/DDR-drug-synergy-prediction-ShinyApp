
# Define server logic required to draw a histogram
server = function(input, output, session) {
  
  #define view names
  shap_view       = "merck_cell_var_synergy_efficacy_21q2_dev_kreis_0.fst"
  inputs_view     = "merck_cell_var_guan_features_21q2_dev_kreis_0.fst"
  prediction_view = "merck_cell_var_guan_predictions_21q2_dev_kreis_0.fst"

  # The model uses a range of features, categorize
  feature_annotation = shiny::reactive({
    create_feature_categories(view = shap_view)
  }) %>%
    shiny::bindCache(shap_view)

  # -- Modules for the SHAP analysis tab
  # These modules create a interactive plot that allows to brose through the
  # different levels of feature categories. There is one for aoc and one for
  # bliss predictions.
  #
  shap_summary_aoc = shiny::reactive({

    # TODO check if there is a translation available for a feature and replace with action name
    synddr::summarize_shap_values(
      feature_annotation = feature_annotation(), 
      shap_view          = shap_view,
      input_view        = inputs_view,
      active_x           = selected_aoc_category(),
      readout            = "aoc")
  }) %>%
    shiny::bindCache(shap_view,
                     selected_aoc_category())

  selected_aoc_category = nested_plot_server(
    id = "aoc", summary = shap_summary_aoc,
    top_level_name = "All Feature Categories",
    ignore_click=c("RDK", "FP2", "Morgan", "MACCS", "FP4", "FP3"))


  shap_summary_bliss = shiny::reactive({

    # TODO check if there is a translation available for a feature and replace with action name
    synddr::summarize_shap_values(
      feature_annotation = feature_annotation(),
      shap_view          = shap_view,
      input_view        = inputs_view,
      active_x           = selected_bliss_category(),
      readout            = "bliss")
  }) %>%
    shiny::bindCache(shap_view,
                     selected_bliss_category())

  selected_bliss_category = nested_plot_server(
    id = "bliss", summary = shap_summary_bliss,
    top_level_name = "All Feature Categories",
    ignore_click=c("RDK", "FP2", "Morgan", "MACCS", "FP4", "FP3"))


  input_summary = shiny::reactive({

    # TODO check if there is a translation available for a feature and replace with action name
    synddr::summarize_input_values(
      feature_annotation = feature_annotation(),
      view               = inputs_view,
      active_x           = selected_input_category())
  }) %>%
    shiny::bindCache(inputs_view,
                     selected_input_category())

  selected_input_category = nested_plot_server(
    id = "input", summary = input_summary,
    top_level_name = "All Feature Categories",
    ignore_click=c("RDK", "FP2", "Morgan", "MACCS", "FP4", "FP3",
                   "Treatment_1_name", "Treatment_2_name", "Treatment_1_moa",
                   "Treatment_2_moa"))


  prediction_server(id      = "aoc_bliss_prediction",
                    view    = prediction_view)

  score_data_server("score_data")
}
