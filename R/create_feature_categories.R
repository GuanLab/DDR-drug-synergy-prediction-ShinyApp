#' @title 
#' Map feature names to categories as used in the publication
#' 
#' @param view A character name of the fiew listing all available features
#' 
create_feature_categories = function(view) {

  # load unique available feature_categories
  feature_categories = load_fst(fst_name = view, cols = "feature", 
                                distinct = TRUE)
  
  # define top level categories
  feature_categories = feature_categories %>%
    dplyr::mutate(category = ifelse(grepl("moa", feature), "Mode-of-action",
                                    ifelse(grepl("gset", feature), "Geneset Annotations",
                                           ifelse(grepl("(_MACCS_|_Morgan_|_RDK_|_FP2_|_FP3_|_FP4_)", feature), "Chemical Structure",
                                                  ifelse(grepl("(_exp$|_cnv$|_lof$|_snv$|_coh_pat$|_lof_pat$|_ddr$)", feature), "Molecular Biomarker",
                                                         ifelse(grepl("Treatment_\\d_name", feature), "Drug Name",
                                                                ifelse(grepl("Geneset", feature), "Geneset",
                                                                       ifelse(grepl("_ave$", feature), "Monotherapy",
                                                                              NA))))))))
  # define molecular subtypes
  feature_categories = feature_categories %>%
    dplyr::mutate(molecular_biomarker = ifelse(grepl("_exp$", feature), "Expression",
                                               ifelse(grepl("ddr$", feature), "DDR",
                                                      ifelse(grepl("coh_pat$", feature), "Coherence",
                                                             ifelse(grepl("cnv$", feature), "CNV",
                                                                    ifelse(grepl("lof_pat$", feature), "LoF",
                                                                           ifelse(grepl("snv$", feature), "SNV",
                                                                                  ifelse(grepl("lof_pat$", feature), "loss-of-function cluster",
                                                                                         ifelse(grepl("lof$", feature), "loss-of-function",
                                                                                                NA)))))))))
  # define fingerprints
  feature_categories = feature_categories %>%
    dplyr::mutate(chemical_structure = ifelse(grepl("MACCS", feature), "MACCS",
                                              ifelse(grepl("FP3", feature), "FP3",
                                                     ifelse(grepl("FP2", feature), "FP2",
                                                            ifelse(grepl("FP4", feature), "FP4",
                                                                   ifelse(grepl("Morgan", feature), "Morgan",
                                                                          ifelse(grepl("RDK", feature), "RDK",
                                                                                 NA)))))))
  # define genes
  feature_categories = feature_categories %>%
    dplyr::mutate(gene = ifelse(grepl("(exp$|cnv$|lof$|snv$|ddr$)", feature), 
                                gsub("_(exp$|cnv$|lof$|snv$|ddr$)", "", feature), 
                                NA))
  feature_categories
}
