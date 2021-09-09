#' #' Features used by the model to predict aoc values
#' #'
#' #' @keywords datasets
#' #' @format A character vector 
#' "aoc_features"
#' 
#' #' Features used by the model to predict bliss values
#' #'
#' #' @keywords datasets
#' #' @format A character vector 
#' "bliss_features"
#' 
#' #' Collection of chemical structure features for all treatments
#' #'
#' #' @keywords datasets
#' #' @format A tibble with a treatment column and the different fingerprints
#' "chemical_structures"
#' 
#' #' Dictionary for drug and affected genes
#' #'
#' #' @keywords datasets
#' #' @format A named list of drug names and hugo symbols
#' "drug2gene"
#' 
#' #' Definition of drugs and their mode of actions
#' #'
#' #' @keywords datasets
#' #' @format A tibble listing drug names (drug_name) and mode of actions (`mode-of-action`)
#' "drug_summary"
#' 
#' #' List of genes used for the surrogate model for aoc prediction
#' #'
#' #' @keywords datasets
#' #' @format A character vector of hugo symbol
#' "genes_aoc"
#' 
#' #' List of genes used for the surrogate model for bliss prediction
#' #'
#' #' @keywords datasets
#' #' @format A character vector of hugo symbols
#' "genes_bliss"
#' 
#' #' Look up table for number of genes hit for each treatment combination
#' #'
#' #' @keywords datasets
#' #' @format A tibble for the treatment combinations and all analysed genesets
#' "geneset_cnt"
#' 
#' #' Table with gene, geneset mapping
#' #'
#' #' @keywords datasets
#' #' @format A tibble with hugo symbols in columns and genesets (1: in geneset, 0: not)
#' "geneset_features"
#' 
#' #' Adjacency matrix with connectivity probabilities
#' #'
#' #' @keywords datasets
#' #' @format A tibble with hugo symbols as column and row names, double indicates connectivity
#' "network"
#' 
#' #' Combination of all drug combination features (except for expressiond ata) used for prediction 
#' #'
#' #' @keywords datasets
#' #' @format A tibble of all required featrures
#' "drug_moa_combn"