#' Build feature set for machine learning model
#' 
#' description
#' 
#' @param mol Molecular signatures from untreated sample
#' @param features drug_name, chemical_structure, molecular (top125), geneset (top125)
#' @param target 'aoc' or bliss'
#' @return 
#' @export
#' 
#' 
<<<<<<< HEAD
build_feature_dataset <- function(data, mol, features, data_path){
=======
##liberary(jsonlite)
build_feature_dataset <- function(data, mol, features, target, data_path){
>>>>>>> 43050be3cb43354dc3034e3a1194536a113e53c9
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
  
  # categorical encoding
  # drug2idx()
  # moa2idx()
  
<<<<<<< HEAD
  build_features <- function(data){
=======
  build_features <- function(row){

>>>>>>> 43050be3cb43354dc3034e3a1194536a113e53c9
    # moa features (2)
    moa2idx(row$)
    # drug name features (2)
    drug2idx()
    # chemical structure feature (xxxx*2)
    all_chemical_structure
    # molecular feature (125)
      #target genes and networks
      target_genes <- unique(c(drug2gene(drug1), drug2gene(drug2)))
      target_genes_exp = paste0(target_genes,'_exp')
      mol[mol[''] == target_genes_exp, ] <-0
      
    #gene set (125) 
      genes <- intersect(colnames(geneset)[-1], target_genes)
      feature_val <- rowSums(geneset[,genes])
      feature_name <- geneset[['geneset']]
<<<<<<< HEAD
    
    # concatennate all features and feature names above
    X, feature_names
  }
  # question: what is the input format for lightGBM prediction? feature names were supposed to be used for SHAP analysis.
  X,feature_names <- build_features(data) 
=======
      
    # X, feature_names
  }
  chemcicals_1 = all_chemical_structure
  colnames(chemcicals_1)[-1] = paste0(colnames(chemcicals_1)[-1], "_1")
  chemcicals_2 = all_chemical_structure
  colnames(chemcicals_2)[-1] = paste0(colnames(chemcicals_2)[-1], "_2")
  
  get_target = 

  apply(data, 1, function(x) {paste0(unlist(drug2gene[[x[[3]]]]), "_exp")})
  
  tmp = data %>%
    dplyr::left_join(chemcicals_1, by=c(".metadata_treatment_1"="treatment")) %>%
    dplyr::left_join(chemcicals_2, by=c(".metadata_treatment_2"="treatment")) %>%
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
      tibble::as.tibble(.)

    res = geneset_counts
    if (length(genes) == 0) {
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

  # X,feature_names <- build_features(data) 
>>>>>>> 43050be3cb43354dc3034e3a1194536a113e53c9
}




