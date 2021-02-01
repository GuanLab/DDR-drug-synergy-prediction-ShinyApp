#' Build feature set for machine learning model
#' @param mol Molecular signatures from untreated sample
#' @param features drug_name, chemical_structure, molecular (top125), geneset (top125)
#' @param target 'aoc' or bliss'
#' @return 
#' @export
#' @importreadr jsonlite
#' 
#' 
##liberary(jsonlite)
build_feature_dataset <- function(data, mol, features, data_path){
  # load feature datasets
  # chemical structure
  all_chemical_structure = read.csv('./data-raw/feature/chemical_structure_features/all_chemical_structure.csv')
  
  # molecular expression data (top 125)
  # ./data_raw/genes_aoc.txt
  # ./data_raw/genes_bliss.txt
  # the union of aoc and bliss top 125 is 182. 
  # check if mol contains all genes and crop the top genes 
  
  
  # gene set (top 125)
  geneset = read.csv('./data-raw/feature/geneset_features/geneset_'+target+'.csv')
  
  # load gene network information from mousenet
  network = jsonlite::read_json('./data-raw/feature/mousenet_features/network.json')  #network[[geneA]][[geneB]]
  
  # load drug-specific target gene information
  drug2gene = jsonlite::read_json('./data-raw/feature/target_gene/drug2gene.json') # example drug2gene[[drugname]] = target gene
  
  # categorical encoding
  drug2idx()
  moa2idx()
  
  build_features <- function(row){
    # moa features (2)
    moa2idx()
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
      
    X, feature_names
  }
 
  X,feature_names <- build_features(data) 
}




