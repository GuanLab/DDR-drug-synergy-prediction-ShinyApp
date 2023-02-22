


test_that("definition of a feature dataset works", {
  mol = readr::read_tsv(system.file("extdata/example_input.tsv",
                                    package = "synddr"),
                        show_col_types = FALSE)
  features = build_feature_dataset(drug_moa_combn, mol) %>%
    dplyr::rename(Treatment_1_moa=.metadata_moa_1,
                  Treatment_2_moa=.metadata_moa_2,
                  Treatment_1_name=.metadata_treatment_1,
                  Treatment_2_name=.metadata_treatment_2)

  bliss_df = features %>%
    tidyr::unnest(bliss_exp)

  testthat::expect_true(all(bliss_features %in% colnames(bliss_df)))

  aoc_df = features %>%
    tidyr::unnest(aoc_exp)

  testthat::expect_true(all(aoc_features %in% colnames(aoc_df)))
})
