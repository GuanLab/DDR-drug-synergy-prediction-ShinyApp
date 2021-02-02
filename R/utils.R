
# create combination of all drugs and MoAs
.combn_pair <- function(x, col) {
  combn(x[[col]], m=2, simplify=F) %>% 
    unlist %>% 
    matrix(., ncol=2, byrow=T) %>%
    data.frame()
}

.access_exp_gset_info <- function(x, y, target, mol){
  
  # identify target genes of drug 1 and 2
  x_targets = unlist(drug2gene[[x]])
  y_targets = unlist(drug2gene[[y]])
  
  # get most important gene names (surogate model)
  gene_sub = gsub("_exp", "", get(paste0("genes_", target)))
  
  # only use the intersection
  genes = intersect(gene_sub, c(x_targets, y_targets))
  
  # count the number of genes that are present in the signatures
  res = get(paste0("geneset_", target))[,c("geneset", genes)] %>% 
    tibble::column_to_rownames("geneset") %>% 
    rowSums %>% 
    t() %>%
    tibble::as_tibble(.)
  
  # if any genes of the surogate model are part of the target genes set
  # expression for these genes to 0 and combine
  # TODO: genes available in mol?
  if (length(genes) > 0) {
    genes = paste0(genes, "_exp")
    mol_ = mol[, get(paste0("genes_", target))]
    mol_[, genes] = 0
    res = cbind(res, mol_)
  }
  
  res
}
