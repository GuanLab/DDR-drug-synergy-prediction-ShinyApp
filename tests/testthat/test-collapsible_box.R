test_that("collapsible box works", {

  shiny::testServer(
    collapsible_box_server,
    args = list(collapsible_title = "Testbox"), {

      testthat::expect_false(open())

      session$setInputs(collapse = 1)

      testthat::expect_true(open())


      session$setInputs(collapse = 1)

      testthat::expect_false(open())


  })


})
