## code to prepare `drug_summary` dataset goes here

drug_summary = readr::read_csv(here::here('inst/extdata/feature/QC/all_drugs_summary.csv'))

usethis::use_data(drug_summary, overwrite = TRUE)
