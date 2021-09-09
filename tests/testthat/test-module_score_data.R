test_that("scoring user uploaded data works", {
  shiny::testServer(
    score_data_server, {
      session$setInputs(file1 = list(datapath = system.file("extdata/example_input.tsv",
                                            package = "guanlabddrdrugcombination")))

      df = readr::read_tsv(system.file("extdata/example_input.tsv",
                                       package = "guanlabddrdrugcombination"),
                           show_col_types = FALSE)

      testthat::expect_equal(valid_df(), df)

      pred = prediction()

      testthat::expect_type(pred, "list")
      testthat::expect_equal(names(pred), "prediction")
      
      efficacy_df = efficacy_stat()
      
      efficacy_df %>%
        .$estimate %>%
        unlist() %>%
        unname() %>%
        testthat::expect_equal(.01535, tolerance = .001)
      
      efficacy_df %>%
        .$p.value %>%
        unlist() %>%
        unname() %>%
        testthat::expect_equal(0.3155, tolerance = .001)
      
      session$setInputs(aoc_bliss_sample = "model 1")
      session$setInputs(nextBtn = 1)
      session$setInputs(prvBtn = 1)
      
      session$setInputs(aoc_bliss_selected = 10)
      session$setInputs(mytable_rows_selected = 8)
    })
})
