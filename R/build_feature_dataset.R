#' Build feature set for machine learning model
#' 
#' description
#' 
#' @param data
#' @param mol Molecular signatures from untreated sample
#' @param target 'aoc' or bliss'
#' 
#' @return 
#' @export
#' 
build_feature_dataset <- function(data, mol, target) {
  
  # load feature datasets
  
  # load gene network information from mousenet
  # TODO: ?
  network = synergyddrapp::network
  
  # add checmical structure features for treatment 1 and treeatment 2
  data <- data %>%
    dplyr::left_join(chemical_structures, 
                     by=c(".metadata_treatment_1"="treatment")) %>%
    dplyr::left_join(chemical_structures, 
                     by=c(".metadata_treatment_2"="treatment"),
                     suffix=c("_1", "_2"))
  
  # add information for target genes and genesets
  data <- data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(target_genes_1 = purrr::map2(.x     = .metadata_treatment_1, 
                                               .y     = .metadata_treatment_2, 
                                               .f     = .access_exp_gset_info,
                                               target = target,
                                               mol    = mol)) %>%
    tidyr::unnest(target_genes_1) 
  
  # recode drug and moa to integers
  data %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(".meta"),
                                ~ as.numeric(as.factor(.x))))
}




