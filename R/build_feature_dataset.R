#' Build feature set for machine learning model
#'
#' description
#'
#' @param data combination of MoA and treatments
#' @param mol tibble Molecular signatures from untreated sample
#'
build_feature_dataset <- function(data, mol) {

  # subset molecular data to surrogate model
  mol_aoc = mol[, genes_aoc] %>%
    dplyr::rename_with(~paste0(.x, "_exp"))
  mol_bliss = mol[, genes_bliss] %>%
    dplyr::rename_with(~paste0(.x, "_exp"))

  adj_mol_aoc = .adjust_by_connectivity(mol_ = mol_aoc, target = "aoc") %>%
    .knockdown_targets(target =  "aoc") %>%
    tidyr::nest(aoc_exp=-dplyr::starts_with("."))

  adj_mol_bliss = .adjust_by_connectivity(mol_bliss, "bliss") %>%
    .knockdown_targets(target = "bliss") %>%
    tidyr::nest(bliss_exp=-dplyr::starts_with("."))

  # encode drug and moa to integers
  id_cols = c(".metadata_moa_1", ".metadata_moa_2", ".metadata_treatment_1",
              ".metadata_treatment_2")

  treatment_map = unique(c(data$.metadata_treatment_1,
                           data$.metadata_treatment_2))
  treatment_map = treatment_map[order(treatment_map)]
  val = 1:length(treatment_map)
  names(val) = treatment_map


  moa_map = unique(c(as.character(data$.metadata_moa_1),
                     as.character(data$.metadata_moa_2)))
  moa_map = moa_map[order(moa_map)]
  moa = 1:length(moa_map)
  names(moa) = moa_map

  # combine all features in same order as for training
  adj_mol_aoc %>%
    dplyr::left_join(data, by = id_cols) %>%
    dplyr::left_join(adj_mol_bliss, by=id_cols) %>%
    dplyr::left_join(geneset_cnt, by = c(".metadata_treatment_1",
                                                ".metadata_treatment_2")) %>%
    dplyr::select(-targets_aoc, -targets_bliss, -network_genes_aoc,
                  -network_genes_bliss, -connectivity_aoc, -connectivity_bliss) %>%
    dplyr::mutate(.metadata_treatment_1 = val[.metadata_treatment_1],
                  .metadata_treatment_2 = val[.metadata_treatment_2],
                  .metadata_moa_1 = moa[.metadata_moa_1],
                  .metadata_moa_2 = moa[.metadata_moa_2]) %>%
    dplyr::select(dplyr::ends_with(c("moa_1", "moa_2")),
                  dplyr::ends_with(c("treatment_1", "treatment_2")),
                  dplyr::starts_with(c("Treatment_1", "Treatment_2")),
                  # dplyr::contains(c("MACCS", "Morgan", "FP2", "FP3", "FP4",
                  #                   "RDK")),
                  dplyr::ends_with("_exp"),
                  dplyr::contains("Geneset"))

  # TODO: check feature order: moa idx, d chemical structure (T1 MACC Morgan,
  # FP2, FP3, FP4, RDK), expressiondata (alpahbetisch),  geneset
}



.adjust_by_connectivity = function(mol_, target) {

  # find available genes
  target_genes = intersect(paste0(get(paste0("genes_", target)), "_exp"),
                           colnames(mol_))

  # find network genes for target genes and combine in a matrix with connectivity
  # scores and NAs
  conn = drug_moa_combn[[glue::glue("connectivity_{target}")]] %>%
    plyr::rbind.fill.matrix()

  # set NAs to 1 (no adjustment of expression values)
  conn[is.na(conn)] = 1
  colnames(conn) = paste0(colnames(conn), "_exp")

  # append all genes that are not related to target genes (no adjustment of
  # expression values)
  missing_genes = setdiff(colnames(mol_), colnames(conn))

  m_ge = matrix(1, nrow=nrow(conn), ncol=length(missing_genes),
                dimnames = list(NULL, missing_genes))

  conn = cbind(conn, m_ge)
  conn = conn[, colnames(mol_)]

  # rowwise multiplication of molecular data and connectivity scores
  mol_ %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::nest(data = -sample) %>%
    dplyr::mutate(data = purrr::map(data,
                                    ~ dplyr::as_tibble(sweep(conn, MARGIN=2,
                                                             unlist(.x), `*`))))
}

.knockdown_targets = function(mol, target) {

  ge = drug_moa_combn[[paste0("targets_", target)]] %>%
    purrr::map(~ matrix(rep(0, length(.x)), nrow=1,
                        dimnames=list(NULL, .x))) %>%
    plyr::rbind.fill.matrix()

  ge[is.na(ge)] = 1
  colnames(ge) = paste0(colnames(ge), "_exp")

  # add genes that aren't targets of any treatment
  missing_genes = setdiff(colnames(mol), colnames(ge))

  m_ge = matrix(1, nrow=nrow(ge), ncol=length(missing_genes),
                dimnames = list(NULL, missing_genes))

  ge = cbind(ge, m_ge)
  ge = ge[, colnames(mol)]

  # set expression to zero and annotation moa/treatment combinations
  mol %>%
    dplyr::mutate(data = purrr::map(data, ~ .x * ge))


  # create multiplication matrix for each sample to set target gene expression
  # to zero
  combination_anno = drug_moa_combn %>%
    dplyr::select(dplyr::starts_with("."))
  mol %>%
    dplyr::mutate(data = purrr::map(data, ~ combination_anno  %>% dplyr::bind_cols(.x))) %>%
    tidyr::unnest(data)
}
