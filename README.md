
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)[![Last-changedate](https://img.shields.io/badge/last%20change-2022--03--18-yellowgreen.svg)](/commits/master)[![minimal
R
version](https://img.shields.io/badge/R%3E%3D-3.6.0-6666ff.svg)](https://cran.r-project.org/)[![packageversion](https://img.shields.io/badge/Package%20version-0.0.2-orange.svg?style=flat-square)](commits/master)

<!-- badges: start -->

[![R-CMD-check](https://github.com/GuanLab/DDR-drug-synergy-prediction-ShinyApp/workflows/R-CMD-check/badge.svg)](https://github.com/GuanLab/DDR-drug-synergy-prediction-ShinyApp/actions)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# guanlabddrdrugcombination

Welcome to our shiny app for DDR drug synergy and efficacy prediction
repository. This repository implements a shiny app that provides a
general overview of all features used for the predictions, the
predictions on out 67 cell lines and summary plots for SHAP analysis on
out data. Additionally, you can upload molecular data and analyse
predicted synergy and efficacy scores of our model

## Requirements

-   R 4.1.1
-   git2r
-   renv

## Installation

To run the app, you will need to install the devtools, git2r and renv
package. The following code will clone and restore the R package
versions that were tested by us.

``` r
install.packages(c('renv', 'git2r'))

# Clone repository
git2r::clone("https://github.com/GuanLab/DDR-drug-synergy-prediction-ShinyApp.git", 
             branch = "renv_support", local_path="guanlabddrdrugcombination")

# Initialize and restore packages using renv
renv::activate(project = "guanlabddrdrugcombination")
renv::restore(project = "guanlabddrdrugcombination", prompt = FALSE)

shiny::runApp('guanlabddrdrugcombination/inst/app/')
```

After a restart of your R session, you can start the app using:

``` r
library(guanlabddrdrugcombination)

launch_synergy_ddr_app()
```
