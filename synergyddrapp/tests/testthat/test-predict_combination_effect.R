test_that("output for optimal combination prediction is equal", {
  
  # access example file
  ex_df <- synergyddrapp::example_input

  # define path to file (required for python access)
  data_path = system.file(package="synergyddrapp")

  # run python implementation
  python_impl <- predict_optimal_drug_combination(mol_df=ex_df, data_path=data_path)

  # run r implementation
  r_impl <- synergyddrapp::predict_optimal_drug_combination(mol_df=ex_df, data_path=data_path)

  # compare returned data frames
  attr(python_impl, 'pandas.index') = NULL
  testthat::expect_equal(r_impl, python_impl)
})