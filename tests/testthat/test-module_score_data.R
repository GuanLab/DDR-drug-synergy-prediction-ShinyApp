test_that("scoring user uploaded data works", {
  shiny::testServer(
    score_data_server, {
      session$setInputs(file1 = list(datapath = system.file("extdata/example_input.tsv",
                                            package = "synddr")))

      df = readr::read_tsv(system.file("extdata/example_input.tsv",
                                       package = "synddr"),
                           show_col_types = FALSE)

      valid_input = shiny::isolate(valid_df())
      testthat::expect_equal(valid_input, df)

      pred = shiny::isolate(prediction())

      testthat::expect_type(pred, "list")
      testthat::expect_equal(names(pred), "prediction")
      
      efficacy_df = shiny::isolate(efficacy_stat())
      testthat::expect_gte(nrow(efficacy_df), 1)
      
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
      
      session$setInputs(aoc_bliss_sample = efficacy_df$sample[[1]])
      session$setInputs(nextBtn = 1)
      session$setInputs(prvBtn = 1)
      
      session$setInputs(aoc_bliss_selected = 10)
      session$setInputs(mytable_rows_selected = 8)
    })
})
