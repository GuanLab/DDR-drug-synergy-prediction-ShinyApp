## code to prepare `drug_summary` dataset goes here

drug_summary = readr::read_csv(here::here('data-raw/feature/QC/all_drugs_summary.csv'))

usethis::use_data(drug_summary, overwrite = TRUE)
