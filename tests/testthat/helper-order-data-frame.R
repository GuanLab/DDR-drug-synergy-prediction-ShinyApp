order_data_frame = function(df, columns = names(df)) {
  ordered = df[do.call(order, unname(df[columns])), , drop = FALSE]
  rownames(ordered) = NULL
  ordered
}

expect_equal_data_frame = function(actual, expected, columns = intersect(names(actual),
                                                                         names(expected))) {
  testthat::expect_equal(order_data_frame(actual, columns),
                         order_data_frame(expected, columns),
                         ignore_attr = TRUE)
}
