
options(fst_folder = "src/fst")

test_that("module for predictions works", {
  
  prediction_view = "merck_cell_var_guan_predictions_21q2_dev_kreis_0.fst"
  
  shiny::testServer(
    prediction_server, 
    args = list(id = "aoc_bliss_prediction",
                view=prediction_view), {
                  
                  session$setInputs(aoc_bliss_sample = "ONCOLEAD_CELL_PC3")
                  
                  pred = predictions()
                  
                  efficacy_df = efficacy_stat()
                  
                  efficacy_df %>%
                    .$estimate %>%
                    unlist() %>%
                    unname() %>%
                    .[[1]] %>%
                    testthat::expect_equal(0.5207632, tolerance = .001)
                  
                  efficacy_df %>%
                    .$p.value %>%
                    unlist() %>%
                    unname() %>%
                    .[[1]] %>%
                    testthat::expect_equal(2.768813e-41, tolerance = .001)
                  
                  session$setInputs(nextBtn = 1)
                  session$setInputs(prvBtn = 1)
                  
                  session$setInputs(aoc_bliss_selected = 10)
                  session$setInputs(mytable_rows_selected = 8)
      
    })
})
