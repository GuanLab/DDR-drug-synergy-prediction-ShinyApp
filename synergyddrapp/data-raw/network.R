## code to prepare `network` dataset goes here

network = jsonlite::read_json(here::here('inst/extdata/feature/mousenet_features/network.json'))

usethis::use_data(network, overwrite = TRUE)
