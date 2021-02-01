## code to prepare `chemical_structures` dataset goes here


genes_aoc = readr::read_tsv(here::here('inst/extdata/feature/molecular/genes_aoc.txt'), col_names=FALSE)[[1]]

usethis::use_data(genes_aoc, overwrite = TRUE)
