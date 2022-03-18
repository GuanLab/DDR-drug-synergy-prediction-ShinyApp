#'
#'
#' Author: julian.kreis@external.merkcgroup.com
#' Date: Tue Feb  1 17:00:33 2021
#' --------------
library(ggplot2)
library(dplyr)
library(DT)
library(ggiraph)

# Set caching options
shiny::shinyOptions(cache = cachem::cache_disk("./myapp-cache"))

shiny::shinyApp(server=guanlabddrdrugcombination:::initialize_server(),
                ui=guanlabddrdrugcombination:::initialize_ui())
