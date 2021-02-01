## code to prepare `example_input` dataset goes here

example_input = readr::read_tsv(here::here('inst/extdata/input/example_input.tsv'))

usethis::use_data(example_input, overwrite = TRUE)
