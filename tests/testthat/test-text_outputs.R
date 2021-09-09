library(shiny)

test_that("Text output functions work", {

  testthat::expect_type(abstract(), "character")
  testthat::expect_type(app_news_list(), "list")
  testthat::expect_type(app_welcome_paragraph(), "character")
})
