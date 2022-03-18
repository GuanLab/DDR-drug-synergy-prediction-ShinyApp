#' Save complete app requirements in compact format
#'
#' Zip all required data for the app. For the app usage no other data than the
#' zip file is needed
#'
#' Author: julian.kreis@external.merkcgroup.com
#' Date: Fri Mar  18 08:48:11 2022
#' --------------

my_wd      = getwd()
dest_path  = here::here('inst/app/')
setwd(dest_path)

files = list.files(dest_path)
files = files[files != "tests"]

zip(zipfile = here::here("inst/extdata/app_data.zip"), 
    files = files)

setwd(my_wd)
