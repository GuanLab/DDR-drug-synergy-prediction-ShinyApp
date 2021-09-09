
#' Create combination of a vector of size m
#' 
#' @param x A data frame
#' @param col A character column name of the data frame x
#' @param m An integer for the number of elements
#' @return A data frame with all combinations with m columns
.combn_pair <- function(x, col, m=2) {
  utils::combn(x[[col]], m=m, simplify=F) %>%
    unlist() %>%
    matrix(., ncol=2, byrow=T) %>%
    data.frame()
}