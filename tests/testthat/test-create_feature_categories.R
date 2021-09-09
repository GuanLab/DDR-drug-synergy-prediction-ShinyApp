library(dplyr)

options(fst_folder = "src/fst")

test_that("defenition of feature categories works", {
  
  view = "merck_cell_var_guan_features_21q2_dev_kreis_0.fst"
  
  feature_annotations = create_feature_categories(view)
  
  testthat::expect_setequal(colnames(feature_annotations),
                            c("feature", "category", "molecular_biomarker",
                              "chemical_structure", "gene"))
  
  testthat::expect_setequal(unique(feature_annotations$chemical_structure),
                            c(NA, "FP2", "FP3", "FP4", "MACCS", "Morgan", "RDK"))
  
  testthat::expect_true(all(genes_aoc %in% unique(feature_annotations$gene)))
  # TODO: not any bliss gene in contribution view?
  # testthat::expect_true(all(genes_bliss %in% unique(feature_annotations$gene)))
  
})
