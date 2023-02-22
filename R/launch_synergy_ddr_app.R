
#' Launcher for 
#'
#' @param fst_folder path to fst files (only for development required)
#'
#' @return
#' Starts a shiny app
#' @export
#'
#' @examples
#' \dontrun{
#' launch_synergy_ddr_app()
#' }
launch_synergy_ddr_app = function(fst_folder=NULL) {
  
  # set options for fst data access
  if (!is.null(fst_folder))
    options(fst_folder = fst_folder)
  
  app_dir = tempdir()# Set caching options
  
  # Copy relevat files
  file.copy(system.file("extdata/app_data.zip", 
                        package = "synddr"), 
            app_dir, recursive=TRUE)
  
  # Unzip files
  utils::unzip(paste(app_dir, "app_data.zip", sep="/"),
               exdir = app_dir)
  
  # start app
  shiny::runApp(app_dir)
}
