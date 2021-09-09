
library(dplyr)

options(fst_folder = "src/fst")

shap_view = "merck_cell_var_synergy_efficacy_21q2_dev_kreis_0.fst"
inputs_view     = "merck_cell_var_guan_features_21q2_dev_kreis_0.fst"

shap_summary_aoc = summarize_shap_values(
  feature_annotation = create_feature_categories(view = shap_view),
  shap_view          = shap_view,
  active_x           = NULL,
  readout            = "aoc",
  input_view         = inputs_view)

test_that("Top level summary works", {

  shiny::testServer(
    nested_plot_server,
    args = list(id = "aoc",
                summary = shiny::reactiveVal(),
                top_level_name = "All Feature Categories"), {
                  
                  summary(shap_summary_aoc)
                  
                  testthat::expect_null(clicked_category())
                  testthat::expect_null(last_last_click())
                  testthat::expect_null(last_clicked_category())

                  testthat::expect_equal(level(), 0)
                  testthat::expect_equal(last_click(),
                                         top_level_name)

                  session$setInputs(categories_plt_selected = "Molecular Biomarker")


                  testthat::expect_null(last_last_click())
                  testthat::expect_null(last_clicked_category())
                  
                  testthat::expect_equal(clicked_category(), "Molecular Biomarker")

                  testthat::expect_equal(level(), 1)
                  testthat::expect_equal(last_click(),
                                         "Molecular Biomarker")
                  
                  shap_summary_aoc = summarize_shap_values(
                    feature_annotation = create_feature_categories(view = shap_view),
                    shap_view          = shap_view,
                    active_x           = "Molecular Biomarker",
                    readout            = "aoc",
                    input_view         = inputs_view)
                  
                  summary(shap_summary_aoc)
                  
                  session$setInputs(categories_plt_selected = "Expression")
                  
                  testthat::expect_equal(last_click(), "Expression")
                  testthat::expect_equal(last_x(), "x")
                  testthat::expect_true(dplyr::all_equal(last_summary(), summary()$df))
                  
                  shap_summary_aoc = summarize_shap_values(
                    feature_annotation = create_feature_categories(view = shap_view),
                    shap_view          = shap_view,
                    active_x           = "Expression",
                    readout            = "aoc",
                    input_view         = inputs_view)
                  
                  summary(shap_summary_aoc)
                  
                  session$setInputs(categories_plt_selected = "ADORA2A_exp")
                  
                  testthat::expect_equal(last_click(), "Expression")
                  testthat::expect_equal(last_x(), "x")
                  
                  session$setInputs(prev = 1)
                  testthat::expect_equal(last_click(), "Molecular Biomarker")
                  
                  
                  session$setInputs(prev = 2)
                  testthat::expect_equal(last_click(), top_level_name)
                })

})
