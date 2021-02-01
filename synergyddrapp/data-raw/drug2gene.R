## code to prepare `network` dataset goes here

drug2gene = jsonlite::read_json(here::here('inst/extdata/feature/target_gene/drug2gene.json'))

usethis::use_data(drug2gene, overwrite = TRUE)
