context('combination effect prediction')

# access example file
ex_df <- synergyddrapp::example_input

reticulate::source_python(here::here("inst/python/main_prediction.py"))

test_that("output of R and python is equal", {
  
  python_impl <- predict_optimal_drug_combination(mol_df=ex_df)
  r_impl      <- synergyddrapp::predict_optimal_drug_combination(mol_df=ex_df)

  cols = colnames(python_impl)[1:4]

  python_impl = python_impl %>%
    dplyr::arrange(!!!rlang::syms(cols))

  r_impl = r_impl %>%
    dplyr::arrange(!!!rlang::syms(cols))

})
