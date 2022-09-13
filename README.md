
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)[![Last-changedate](https://img.shields.io/badge/last%20change-2021--09--09-yellowgreen.svg)](/commits/master)[![minimal
R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)[![packageversion](https://img.shields.io/badge/Package%20version-0.0.0.9000-orange.svg?style=flat-square)](commits/master)

[![R-CMD-check](https://github.com/GuanLab/DDR-drug-synergy-prediction-ShinyApp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GuanLab/DDR-drug-synergy-prediction-ShinyApp/actions/workflows/R-CMD-check.yaml)

# DeepDDR

Welcome to our shiny app for DDR drug synergy and efficacy prediction
repository. This repository implements a shiny app that provides a
general overview of all features used for the predictions, the
predictions on our 62 cell lines and summary plots for SHAP analysis on
our data. Additionally, a user has the opportunity to upload molecular data and analyse
predicted synergy and efficacy scores of our model

<img width="1438" alt="App screenshot for predicted efficacy and synergy scores as described in publication" src="https://user-images.githubusercontent.com/25581180/132656337-80ced6ac-e184-4188-aefb-42006f902119.png">

For scoring new experimental data, the user only needs to upload the TPM expression values for 78 genes:

<img width="1437" alt="Screenshot for predicting efficacy and synergy scores" src="https://user-images.githubusercontent.com/25581180/132656155-5f504c6f-db81-4a64-be9f-897b1573e6bd.png">


## Requirements

This application was tested using R 4.1.1.

## Installation

To run the app, you will need to install the devtools package.

``` r
install.packages('devtools')
```

Next you can install the app using:

``` r
devtools::install_github("GuanLab/DDR-drug-synergy-prediction-ShinyApp")
```

After a restart of your R session, you can start the app using:

``` r
library(guanlabddrdrugcombination)

launch_synergy_ddr_app()
```
