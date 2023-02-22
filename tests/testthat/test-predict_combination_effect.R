test_that("predict combination effect works", {
  
  mol = readr::read_tsv(system.file("extdata", "example_input.tsv",
                                    package = "synddr"),
                        show_col_types = FALSE)
  
  features <- build_feature_dataset(data=drug_moa_combn,
                                    mol=mol)
  
  aoc_X = features %>%
    dplyr::select(-dplyr::contains("bliss")) %>%
    tidyr::unnest(dplyr::contains("aoc"))
  
  sample_features = aoc_X %>%
    dplyr::select(sample, dplyr::starts_with("."))
  
  aoc_X = aoc_X %>%
    dplyr::select(-sample) %>%
    dplyr::rename(Treatment_1_moa=.metadata_moa_1,
                  Treatment_2_moa=.metadata_moa_2,
                  Treatment_1_name=.metadata_treatment_1,
                  Treatment_2_name=.metadata_treatment_2) %>%
    dplyr::select(!!!aoc_features)
  
  # predict efficacy: Not enough memory available for github actions on windows
  # aoc_Y = predict_combination_effect(X           = aoc_X,
  #                                    pred_target = 'aoc',
  #                                    path        = system.file("models",
  #                                                              package = "synddr"),
  #                                    threads     = 1,
  #                                    predcontrib = TRUE)
})
