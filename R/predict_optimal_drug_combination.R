#' Predict best drug combinations based on patient's molecular readouts
#'
#' @param mol_df A data.frame with molecular readouts from patient
#' @param predcontrib A logical if per-feature contributions should be calculated
#' @param model_path A character to the path containing model files
#' (name: '{pred_target}_{cnt}.model.txt')
#' @param threads An integer for the number of threads used for parallelization
#' @param chunk_size An integer defining the chunck size that is 
#'
#' @return best drug combinations based on prediction
#'
predict_optimal_drug_combination = function(mol_df, predcontrib=FALSE,
                                            threads=1, model_path=NULL,
                                            chunk_size = 50) {
  shiny::withProgress({

    mol_df = mol_df %>%
      tibble::column_to_rownames("sample")

    # chunk the data into subsets of size chunk_size
    subs    = seq(0, nrow(mol_df)-1, by=chunk_size)
    sub_dfs = purrr::map(subs, ~ mol_df[(.x + 1):(min(nrow(mol_df), .x + chunk_size)), ])

    shiny::incProgress(.1, message = "Score data",
                       detail = glue::glue("Process {length(sub_dfs)} chunks"))

    chunck_step_size = .9 / length(sub_dfs)

    # iterate over chunks
    res = purrr::map(seq_along(sub_dfs), function(idx) {

      mol = sub_dfs[[idx]]

      result = list()

      shiny::incProgress(chunck_step_size / 3,
                         detail =  glue::glue("Chunk {idx/length(sub_dfs)} -",
                                              "Building features"))

      # create a set of features
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

      shiny::incProgress(chunck_step_size / 3,
                         detail =  glue::glue("Chunk {idx} of {length(sub_dfs)} - ",
                                              "Predict efficacy scores"))
      # predict efficacy:
      aoc_Y = predict_combination_effect(X           = aoc_X,
                                         pred_target = 'aoc',
                                         path        = model_path,
                                         threads     = threads,
                                         predcontrib = predcontrib)
      pred_aoc = aoc_Y
      if (predcontrib) {
        pred_aoc = purrr::map(aoc_Y, ~rowSums(.x)) %>%
          Reduce(f = rbind) %>%
          matrix(nrow=length(aoc_Y)) %>%
          colMeans()
        result$shap_aoc = aoc_Y
      }

      # create a set of features
      bliss_X = features %>%
        dplyr::select(-dplyr::contains("aoc")) %>%
        tidyr::unnest(dplyr::contains("bliss")) %>%
        dplyr::select(-sample, -.metadata_treatment_1, -.metadata_treatment_2) %>%
        dplyr::rename(Treatment_1_moa=.metadata_moa_1,
                      Treatment_2_moa=.metadata_moa_2) %>%
        dplyr::select(!!!bliss_features)

      shiny::incProgress(chunck_step_size / 3,
                         detail =  glue::glue("Chunk {idx} of {length(sub_dfs)} - ",
                                              "Predict synergy"))
      bliss_Y = predict_combination_effect(X           = bliss_X,
                                           pred_target = 'bliss',
                                           path        = model_path,
                                           threads     = threads,
                                           predcontrib = predcontrib)
      pred_bliss = bliss_Y
      if (predcontrib) {
        pred_bliss = purrr::map(bliss_Y, ~ rowSums(.x)) %>%
          Reduce(f = rbind) %>%
          matrix(nrow=length(bliss_Y)) %>%
          colMeans()

        result$shap_bliss = bliss_Y
      }
      predictions = sample_features %>%
        dplyr::mutate(AOC    = pred_aoc,
                      Bliss  = pred_bliss)

      treatment_map = unique(c(drug_moa_combn$.metadata_treatment_1,
                               drug_moa_combn$.metadata_treatment_1))
      treatment_map = treatment_map[order(treatment_map)]
      val = 1:length(treatment_map)
      names(val) = treatment_map



      moa_map = unique(c(drug_moa_combn$.metadata_moa_1,
                         drug_moa_combn$.metadata_moa_2))
      moa_map = moa_map[order(moa_map)]
      moa = 1:length(moa_map)
      names(moa) = moa_map

      shiny::incProgress(, detail = glue::glue("Chunk {idx} of {length(sub_dfs)} - ",
                                                                   "Predict efficacy"))
      # append predictions to treatment/moa combinations
      result$prediction <- drug_moa_combn %>%
        dplyr::select(dplyr::starts_with(".")) %>%
        dplyr::transmute(dplyr::across(dplyr::starts_with("."), ~ .x,
                                       .names = "{gsub('\\\\.', '', col)}")) %>%
        dplyr::mutate(.metadata_treatment_1 = val[metadata_treatment_1],
                      .metadata_treatment_2 = val[metadata_treatment_2],
                      .metadata_moa_1 = moa[metadata_moa_1],
                      .metadata_moa_2 = moa[metadata_moa_2]) %>%
        dplyr::left_join(predictions, by = c(".metadata_moa_1",
                                             ".metadata_moa_2",
                                             ".metadata_treatment_1",
                                             ".metadata_treatment_2")) %>%
        dplyr::select(-dplyr::starts_with(".")) %>%
        dplyr::rename(`MoA 1`=metadata_moa_1, `MoA 2`=metadata_moa_2,
                      `Treatment 1`=metadata_treatment_1,
                      `Treatment 2`=metadata_treatment_2) %>%
        dplyr::filter(!is.na(AOC), !is.na(Bliss)) # TODO: Fix categorisation using python implementation

      result
    })

    # merge chunks
    items = names(res[[1]])
    dfs = purrr::map_lgl(res[[1]], ~ "tbl" %in% class(.x))

    res1 =purrr::map(items[dfs], function(item) purrr::map_df(res, ~ .x[[item]]))
    names(res1) = items[dfs]

    lists = purrr::map_lgl(res[[1]], ~ "list" %in% class(.x))

    res2 =purrr::map(items[lists], function(item) {
      purrr::map(seq_along(res[[1]][[item]]), function(fold) {
        purrr::map_df(res, ~ .x[[item]][[fold]])
      })
    })
    names(res2) = items[lists]

    return(append(res1, res2))
  })
}
