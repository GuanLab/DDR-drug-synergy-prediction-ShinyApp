
library(dplyr)

options(fst_folder = "src/fst")

shap_view = "merck_cell_var_synergy_efficacy_21q2_dev_kreis_0.fst"
input_view = "merck_cell_var_guan_features_21q2_dev_kreis_0.fst"

feature_annotation = create_feature_categories(view = shap_view)

test_that("Level 1 categories are summarized correctely", {

  molecular_biomarker_features = feature_annotation %>%
    dplyr::pull(feature)

  all_aoc = summarize_shap_values(feature_annotation, shap_view = shap_view, 
                                  input_view = input_view, NULL, "aoc")

  all_bliss = summarize_shap_values(feature_annotation, shap_view,
                                    input_view = input_view,  NULL, "bliss")

  shaps = load_fst(fst_name=shap_view, 
                   filter = list(feature = molecular_biomarker_features),
                   cols = c(".identifier_sample_name", "feature",
                            paste0("aoc_y_", 1:5))) %>%
    dplyr::left_join(feature_annotation, by="feature") %>%
    dplyr::group_by(.identifier_sample_name, category) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("aoc")),
                        ~ abs(sum(.x, na.rm = TRUE))) %>%
    dplyr::group_by(category) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("aoc")),
                        ~ mean(.x, na.rm = TRUE)) %>%
    tidyr::gather(-category, key="measure", value="shap") %>%
    tidyr::separate(measure, into = c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x=category)

  testthat::expect_true(dplyr::all_equal(all_aoc$df, shaps))

  shaps = load_fst(shap_view, filter = list(feature = molecular_biomarker_features),
                   cols = c(".identifier_sample_name", "feature",
                            paste0("bliss_y_", 1:5))) %>%
    dplyr::left_join(feature_annotation, by="feature") %>%
    dplyr::group_by(.identifier_sample_name, category) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("bliss")),
                        ~ abs(sum(.x, na.rm = TRUE))) %>%
    dplyr::group_by(category) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("bliss")),
                        ~ mean(.x, na.rm = TRUE)) %>%
    tidyr::gather(-category, key="measure", value="shap") %>%
    tidyr::separate(measure, into = c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x=category)

  testthat::expect_true(dplyr::all_equal(all_bliss$df, shaps))
})

test_that("Level 2 categories are summarized correctely", {

  molecular_biomarker_features = feature_annotation %>%
    dplyr::filter(category == "Molecular Biomarker") %>%
    dplyr::pull(feature)

  all_aoc = summarize_shap_values(feature_annotation, shap_view = shap_view, 
                                  input_view = input_view,
                                  "Molecular Biomarker", "aoc")

  all_bliss = summarize_shap_values(feature_annotation, shap_view = shap_view, 
                                    input_view = input_view, 
                                    "Molecular Biomarker", "bliss")

  molecular_biomarker_features = feature_annotation %>%
    dplyr::filter(category == "Molecular Biomarker")
  
  shaps = load_fst(shap_view, filter = list(feature = molecular_biomarker_features$feature),
                   cols = c(".identifier_sample_name", "feature",
                            paste0("aoc_y_", 1:5))) %>%
    dplyr::left_join(feature_annotation, by="feature") %>%
    dplyr::group_by(.identifier_sample_name, molecular_biomarker) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("aoc")),
                        ~ abs(sum(.x, na.rm = TRUE))) %>%
    dplyr::group_by(molecular_biomarker) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("aoc")),
                        ~ mean(.x, na.rm = TRUE)) %>%
    tidyr::gather(-molecular_biomarker, key="measure", value="shap") %>%
    tidyr::separate(measure, into = c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x=molecular_biomarker)

  testthat::expect_true(dplyr::all_equal(all_aoc$df, shaps))

  shaps = load_fst(shap_view, filter = list(feature = molecular_biomarker_features$feature),
                   cols = c(".identifier_sample_name", "feature",
                            paste0("bliss_y_", 1:5))) %>%
    dplyr::left_join(feature_annotation, by="feature") %>%
    dplyr::group_by(.identifier_sample_name, molecular_biomarker) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("bliss")),
                        ~ abs(sum(.x, na.rm = TRUE))) %>%
    dplyr::group_by(molecular_biomarker) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("bliss")),
                        ~ mean(.x, na.rm = TRUE)) %>%
    tidyr::gather(-molecular_biomarker, key="measure", value="shap") %>%
    tidyr::separate(measure, into = c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x = molecular_biomarker)

  testthat::expect_true(dplyr::all_equal(all_bliss$df, shaps))
})

# 
test_that("Level 3 numerical features are summarized correctly", {

  columns = load_fst(fst_name = shap_view, peek=TRUE)

  molecular_biomarker_features = feature_annotation %>%
    dplyr::filter(molecular_biomarker == "Expression") %>%
    dplyr::pull(feature)

  all_aoc = summarize_shap_values(feature_annotation = feature_annotation,
                                  shap_view = shap_view,
                                  input_view = input_view, 
                                  active_x = "Expression",
                                  readout = "aoc")

  all_bliss = summarize_shap_values(feature_annotation, shap_view = shap_view,
                                    input_view = input_view,
                                    "Expression", "bliss")
  cols = columns %>%
    dplyr::select(dplyr::starts_with(c(".iden", ".meta", ".response",
                                       "aoc")),
                  "feature") %>%
    base::colnames()

  shaps = load_fst(shap_view, filter = list(feature = molecular_biomarker_features),
                   cols = cols) %>%
    dplyr::group_by(.identifier_sample_name, feature) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("aoc")),
                        ~ abs(sum(.x, na.rm = TRUE))) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("aoc")),
                        ~ mean(.x, na.rm = TRUE)) %>%
    tidyr::gather(-feature, key="measure", value="shap") %>%
    tidyr::separate(measure, into = c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x = feature)

  testthat::expect_true(dplyr::all_equal(all_aoc$df, shaps))

  cols = columns %>%
    dplyr::select(dplyr::starts_with(c(".iden", ".meta", ".response",
                                       "bliss")),
                  "feature") %>%
    base::colnames()

  shaps = load_fst(shap_view, filter = list(feature = molecular_biomarker_features),
                   cols = cols) %>%
    dplyr::group_by(.identifier_sample_name, feature) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("bliss")),
                        ~ abs(sum(.x, na.rm = TRUE))) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise_at(vars(dplyr::starts_with("bliss")),
                        ~ mean(.x, na.rm = TRUE)) %>%
    tidyr::gather(-feature, key="measure", value="shap") %>%
    tidyr::separate(measure, into = c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x = feature)

  testthat::expect_true(dplyr::all_equal(all_bliss$df, shaps))
})


test_that("Level 3 categorical features are summarized correctely", {

  columns = load_fst(shap_view, peek=TRUE)

  molecular_biomarker_features = feature_annotation %>%
    dplyr::filter(chemical_structure == "FP2") %>%
    dplyr::pull(feature)

  all_aoc = summarize_shap_values(feature_annotation = feature_annotation,
                                  shap_view = shap_view,
                                  input_view = input_view,
                                  active_x = "FP2",
                                  readout = "aoc")

  all_bliss = summarize_shap_values(feature_annotation, shap_view,
                                    input_view = input_view,
                                    "FP2", "bliss")

  cols = columns %>%
    dplyr::select(dplyr::starts_with(c(".iden", ".meta", ".response",
                                       "aoc")),
                  "feature") %>%
    base::colnames()

  shaps = load_fst(shap_view, filter = list(feature = molecular_biomarker_features),
                   cols = cols) %>%
    dplyr::filter(feature != "base_value") %>%
    dplyr::group_by(.identifier_sample_name, feature) %>%
    dplyr::summarize_at(vars(dplyr::starts_with("aoc")),
                        ~ abs(sum(., na.rm=TRUE))) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize_at(vars(dplyr::starts_with("aoc")),
                        ~ mean(., na.rm=TRUE)) %>%
    tidyr::gather(-feature, key="measure", value="shap") %>%
    tidyr::separate(measure, into=c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x = feature)


  testthat::expect_true(dplyr::all_equal(all_aoc$df, shaps))

  cols = columns %>%
    dplyr::select(dplyr::starts_with(c(".iden", ".meta", ".response",
                                       "bliss")),
                  "feature") %>%
    base::colnames()

  shaps = load_fst(shap_view, filter = list(feature = molecular_biomarker_features),
                   cols = cols) %>%
    dplyr::filter(feature != "base_value") %>%
    dplyr::group_by(.identifier_sample_name, feature) %>%
    dplyr::summarize_at(vars(dplyr::starts_with("bliss")),
                        ~ abs(sum(., na.rm=TRUE))) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize_at(vars(dplyr::starts_with("bliss")),
                        ~ mean(., na.rm=TRUE)) %>%
    tidyr::gather(-feature, key="measure", value="shap") %>%
    tidyr::separate(measure, into=c("output", "fold"), sep="_y_") %>%
    dplyr::rename(x = feature)

  testthat::expect_true(dplyr::all_equal(all_bliss$df, shaps))
})
