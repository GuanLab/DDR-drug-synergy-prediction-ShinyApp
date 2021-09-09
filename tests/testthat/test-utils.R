test_that("definition of drug combinations works", {
  
  comb = function(n, x) {factorial(n) / factorial(n-x) / factorial(x)}
  
  res = .combn_pair(drug_summary, "drug_name")
  
  testthat::expect_equal(nrow(res), comb(nrow(drug_summary), 2))
  
})
