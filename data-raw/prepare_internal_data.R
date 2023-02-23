library(dplyr)
library(synddr)


aoc_features = readr::read_tsv(here::here("data-raw/feature/in_aoc_model.tsv"),
                               col_names = FALSE, show_col_types = FALSE) %>%
  unlist() %>%
  unname()

bliss_features = readr::read_tsv(here::here("data-raw/feature/in_bliss_model.tsv"),
                                 col_names = FALSE, show_col_types = FALSE) %>%
  unlist() %>%
  unname()


chemical_structures = readr::read_csv(here::here('data-raw/feature/all_chemical_structure.csv'))

drug2gene = jsonlite::read_json(here::here('data-raw/feature/drug2gene.json'))

drug_summary = readr::read_csv(here::here('data-raw/feature/all_drugs_summary.csv'))

tratment_comb_2 <- synddr:::.combn_pair(drug_summary, "drug_name") %>%
  dplyr::rename(.metadata_treatment_1=X1, .metadata_treatment_2=X2)

moa_comb_2 <- synddr:::.combn_pair(drug_summary, "mode-of-action") %>%
  dplyr::rename(.metadata_moa_1=X1, .metadata_moa_2=X2)

drug_moa_combn <- cbind(moa_comb_2, tratment_comb_2) %>%
  dplyr::as_tibble()

genes_aoc = readr::read_tsv(here::here('data-raw/feature/genes_aoc.txt'),
                            col_names=FALSE, show_col_types = FALSE)[[1]]

genes_aoc = gsub("_exp", "", genes_aoc)


genes_bliss = readr::read_tsv(here::here('data-raw/feature/genes_bliss.txt'), 
                              col_names=FALSE, show_col_types = FALSE)[[1]]

genes_bliss = gsub("_exp", "", genes_bliss)

# add drug target genes
drugs <- drug_moa_combn %>%
  dplyr::rowwise() %>%
  dplyr::mutate(targets = list(c(unlist(drug2gene[[.metadata_treatment_1]]),
                                 unlist(drug2gene[[.metadata_treatment_2]])))) %>%
  dplyr::mutate(targets_aoc = list(intersect(genes_aoc, targets)),
                targets_bliss = list(intersect(genes_bliss, targets))) %>%
  dplyr::select(-targets)

nw = jsonlite::read_json(here::here('data-raw/feature/network.json'))

# select elements with entries
has_interactions = purrr::map_lgl(nw, ~ length(.x) >0)

nw = nw[has_interactions]

net = purrr::map(nw, function(x) {as.data.frame(x)}) %>%
  dplyr::bind_rows()

rownames(net) = names(nw)

network = net

av_genes_aoc = intersect(rownames(network), genes_aoc)
av_genes_bliss = intersect(rownames(network), genes_bliss)


.max_connectivity = function(x, y) {
  ge = network[x, y, drop=F]
  
  rows = apply(ge, 1, function(x) any(!is.na(x)))
  
  res = NULL
  if (any(rows)) {
    ge = ge[rows, , drop=F]
    res = ge %>%
      apply(., 1, max, na.rm=TRUE)
  }
  return(res)
}

# add connectivity scores
drugs = drugs %>%
  dplyr::mutate(network_genes_aoc = list(intersect(av_genes_aoc, targets_aoc)),
                network_genes_bliss = list(intersect(av_genes_bliss, targets_bliss))) %>%
  dplyr::mutate(connectivity_aoc = list(.max_connectivity(av_genes_aoc, network_genes_aoc)),
                connectivity_bliss = list(.max_connectivity(av_genes_bliss, network_genes_bliss))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(connectivity_aoc = purrr::map(connectivity_aoc, function(x) {if(!is.null(x)) {matrix(1-x, nrow=1, dimnames = list(NULL, names(x)))} else {NA}}),
                connectivity_bliss = purrr::map(connectivity_bliss, function(x) {if(!is.null(x)) {matrix(1-x, nrow=1, dimnames = list(NULL, names(x)))} else {NA}}))


# add checmical structure features for treatment 1 and treatment 2
drug_moa_combn <- drugs %>%
  dplyr::ungroup() %>%
  dplyr::left_join(chemical_structures,
                   by=c(".metadata_treatment_1"="treatment")) %>%
  dplyr::left_join(chemical_structures,
                   by=c(".metadata_treatment_2"="treatment"),
                   suffix=c("_1", "_2"))

t1 = drug_moa_combn %>%
  dplyr::select(-dplyr::starts_with(c(".", "target", "network", "connectivity"))) %>%
  dplyr::rename_with(~paste0("Treatment_1_", gsub("_1$", "", .x)),
                     .cols = dplyr::ends_with("_1"))
t2 =drug_moa_combn %>%
  dplyr::select(-dplyr::starts_with(c(".", "target", "network", "connectivity"))) %>%
  dplyr::rename_with(~paste0("Treatment_2_", gsub("_2$", "", .x)),
                     .cols = dplyr::ends_with("_2"))
drug_moa_combn = drug_moa_combn %>%
  dplyr::select(dplyr::starts_with(c(".", "target", "network", "connectivity"))) %>%
  dplyr::bind_cols(t1, t2)

geneset_features = readr::read_csv(here::here('data-raw/feature/geneset_features.csv'))  %>%
  dplyr::rename(geneset=`...1`)

.get_gset_affected = function() {
  drug_moa_combn %>%
    dplyr::select(.metadata_treatment_1, .metadata_treatment_2) %>%
    dplyr::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gset = purrr::map2(.metadata_treatment_1, .metadata_treatment_2,
                                     function(.x, .y) {
                                       append(drug2gene[[.x]], drug2gene[[.y]]) %>%
                                         unlist() %>%
                                         intersect(colnames(geneset_features),
                                                   .) %>%
                                         list()
                                     })) %>%
    dplyr::mutate(gset = purrr::map(gset, .f = function(.x) {
      geneset_features[, c("geneset", .x)] %>%
        tibble::column_to_rownames("geneset") %>%
        rowSums %>%
        t() %>%
        tibble::as_tibble(.)
    })) %>%
    tidyr::unnest(gset) %>%
    tidyr::gather(-dplyr::starts_with("."), key="signature", value="cnt")
  
}

geneset_cnt = .get_gset_affected()

geneset_cnt = geneset_cnt %>%
  dplyr::mutate(signature = paste("Geneset", signature, sep="_"))  %>%
  tidyr::spread(key="signature", value="cnt")


usethis::use_data(aoc_features, bliss_features, chemical_structures,
                  drug2gene, drug_moa_combn, genes_aoc, drug_summary,
                  genes_bliss, geneset_cnt, geneset_features, network,
                  overwrite = TRUE, internal = TRUE)
