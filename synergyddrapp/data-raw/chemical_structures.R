## code to prepare `chemical_structures` dataset goes here


chemical_structures = readr::read_csv(here::here('inst/extdata/feature/chemical_structure_features/all_chemical_structure.csv'))

usethis::use_data(chemical_structures, overwrite = TRUE)
