
library(dplyr)

options(fst_folder = "src/fst")

view = "merck_cell_var_guan_features_21q2_dev_kreis_0.fst"

feature_annotation = create_feature_categories(view = view)

categories = feature_annotation %>%
  dplyr::select(-gene)

# read features from trained lightGBM model
aoc_features = readr::read_tsv("src/in_aoc_model.tsv",
                               show_col_types = FALSE) %>%
  unlist() %>%
  unname()
bliss_features = readr::read_tsv("src/in_bliss_model.tsv",
                                 show_col_types = FALSE) %>%
  unlist() %>%
  unname()
all_features = unique(c(bliss_features, aoc_features))

test_that("Level 1 categories are summarized correctely", {


  level_1 = summarize_input_values(feature_annotation = feature_annotation,
                                   view = view, active_x = NULL)


  # level_1$df %>%
  #   dplyr::distinct(category, Total) %>%
  #   dplyr::pull(Total) %>%
  #   sum() %>%
  #   testthat::expect_equal(length(all_features))

  testthat::expect_equal(level_1$x, "category")
  testthat::expect_equal(level_1$x_title, "All Feature Categories")
  testthat::expect_equal(level_1$y, "x")
  testthat::expect_equal(level_1$y_title, "Number of Features")
  testthat::expect_equal(level_1$fill, "subcategory")
  testthat::expect_equal(level_1$legend_title, "Subcategory")

  # level_1$df %>%
  #   dplyr::filter(category == "Chemical Structure") %>%
  #   dplyr::distinct(category, Total) %>%
  #   dplyr::pull(Total) %>%
  #   testthat::expect_equal(length(grep("(Morgan|RDK|FP2|FP3|FP4|MACCS)", all_features)))

  # level_1$df %>%
  #   dplyr::filter(category == "Molecular Biomarker") %>%
  #   dplyr::distinct(category, Total) %>%
  #   dplyr::pull(Total) %>%
  #   testthat::expect_equal(length(grep("(_exp$)", all_features, value=TRUE)))

  # level_1$df %>%
  #   dplyr::filter(category == "Geneset") %>%
  #   dplyr::distinct(category, Total) %>%
  #   dplyr::pull(Total) %>%
  #   testthat::expect_equal(length(grep("(^Geneset_)", all_features, value=TRUE)))

  # level_1$df %>%
  #   dplyr::filter(category == "Drug Name") %>%
  #   dplyr::distinct(category, Total) %>%
  #   dplyr::pull(Total) %>%
  #   testthat::expect_equal(length(grep("(_name$)", all_features, value=TRUE)))

  # level_1$df %>%
  #   dplyr::filter(category == "Mode-of-action") %>%
  #   dplyr::distinct(category, Total) %>%
  #   dplyr::pull(Total) %>%
  #   testthat::expect_equal(length(grep("(_moa$)", all_features, value=TRUE)))


})

test_that("Level 2 subcategory is summarized correctely", {

  category = "Molecular Biomarker"

  level_2_bm = summarize_input_values(feature_annotation, view,
                                      active_x = category)

  # level_2_bm$df %>%
  #   dplyr::filter(molecular_biomarker == "Expression") %>%
  #   dplyr::distinct(molecular_biomarker, n) %>%
  #   dplyr::pull(n) %>%
  #   testthat::expect_equal(length(grep("(_exp$)", all_features, value=TRUE)))

})

test_that("Level 2 subcategories are summarized correctely", {

  category = "Chemical Structure"

  level_2_cs = summarize_input_values(feature_annotation, view,
                                      active_x = category)

  # purrr::map_lgl(c("RDK", "FP2", "FP3", "FP4", "MACCS"),
  #                function(id) {
  #                  level_2_cs$df %>%
  #                    dplyr::filter(chemical_structure == id) %>%
  #                    dplyr::pull(n) %>%
  #                    `==`(length(grep(glue::glue("(_{id}_)"),
  #                                     all_features, value=TRUE)))
  #                  }) %>%
  #   all() %>%
  #   testthat::expect_true()

})


test_that("Level 2 categorical features are summarized correctely", {

  category = "Geneset"

  level_2_gs = summarize_input_values(feature_annotation, view,
                                      active_x = category)

  testthat::expect_equal(level_2_gs$x, "feature")
  testthat::expect_equal(level_2_gs$x_title, "Geneset")
  testthat::expect_equal(level_2_gs$y, "count")
  testthat::expect_equal(level_2_gs$y_title, "Number of Features")
  testthat::expect_equal(level_2_gs$fill, "x")
  testthat::expect_equal(level_2_gs$legend_title, "Value")

  # level_2_gs$df %>%
  #   dplyr::pull(feature) %>%
  #   as.character() %>%
  #   testthat::expect_setequal(grep("^Geneset", all_features, value = TRUE))
})


test_that("level 3 dumerical features are summarized correctely", {

  category = "Expression"

  level_3_exp = summarize_input_values(feature_annotation, view,
                                      active_x = category)

  testthat::expect_equal(level_3_exp$x, "feature")
  testthat::expect_equal(level_3_exp$x_title, "Expression")

  # level_3_exp$df %>%
  #   dplyr::pull(feature) %>%
  #   testthat::expect_setequal(grep("_exp$", all_features, value = TRUE))

  level_3_exp$df %>%
    colnames() %>%
    testthat::expect_setequal(c("min", "max", "middle", "lower", "upper", "feature"))
})