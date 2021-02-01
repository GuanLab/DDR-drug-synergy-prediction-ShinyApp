#' Build feature set for machine learning model
#' 
#' description
#' 
#' @param mol Molecular signatures from untreated sample
#' @param target 'aoc' or bliss'
#' @return 
#' @export
#' 
#' 
##liberary(jsonlite)
build_feature_dataset <- function(data, mol, target, data_path){
  # load feature datasets
  # chemical structure
  all_chemical_structure = synergyddrapp::chemical_structures
  
  # molecular expression data (top 125)
  # ./data_raw/genes_aoc.txt
  # ./data_raw/genes_bliss.txt
  # the union of aoc and bliss top 125 is 182. 
  # check if mol contains all genes and crop the top genes 
  
  
  # gene set (top 125)
  geneset = get(paste0("geneset_", target))

  # load gene network information from mousenet
  network = synergyddrapp::network  #network[[geneA]][[geneB]]
  
  # load drug-specific target gene information
  drug2gene = synergyddrapp::drug2gene # example drug2gene[[drugname]] = target gene
  
  data %>%
    dplyr::left_join(all_chemical_structure, by=c(".metadata_treatment_1"="treatment")) %>%
    dplyr::left_join(all_chemical_structure, by=c(".metadata_treatment_2"="treatment"),
                     suffix=c("_1", "_2")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(target_genes_1 = purrr::map2(.metadata_treatment_1, .metadata_treatment_2, function(x, y){
    
    x_targets = unlist(drug2gene[[x]])
    y_targets = unlist(drug2gene[[y]])

    gene_sub = gsub("_exp", "", get(paste0("genes_", target)))

    genes = intersect(gene_sub, c(x_targets, y_targets))

    geneset_counts = geneset[,c("geneset", genes)] %>% 
      tibble::column_to_rownames("geneset") %>% 
      rowSums %>% 
      t() %>%
      tibble::as_tibble(.)

    res = geneset_counts
    if (length(genes) > 0) {
      genes = paste0(genes, "_exp")
      mol_ = mol[, get(paste0("genes_", target))]
      mol_[, genes] = 0
      res = cbind(res, mol_)
    }

    res
  })) %>%
  tidyr::unnest(target_genes_1) %>%
    dplyr::mutate(.metadata_moa_1 = as.numeric(as.factor(.metadata_moa_1)),
                  .metadata_moa_2 = as.numeric(as.factor(.metadata_moa_2)),
                  .metadata_treatment_1 = as.numeric(as.factor(.metadata_treatment_1)),
                  .metadata_treatment_2 = as.numeric(as.factor(.metadata_treatment_2))) 

}




