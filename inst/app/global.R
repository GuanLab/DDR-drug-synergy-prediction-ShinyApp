#'
#'
#' Author: julian.kreis@external.merkcgroup.com
#' Date: Tue Feb  1 17:00:33 2021
#' --------------
library(ggplot2)
library(dplyr)
library(DT)
library(ggiraph)
library(synddr)

# set options for fst data access
options(fst_folder = "~/fst")

# Set caching options
shiny::shinyOptions(cache = cachem::cache_disk("./myapp-cache"))
