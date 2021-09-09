
#' Summary function for SHAP values
#'
#' This function summarizes the SHAP values of the test data for a defined
#' set of features
#'
#' @param feature_annotation data.frame defining categories of features (see 
#' [create_feature_categories()])
#' @param shap_view A character name of the fst file that contains SHAP values
#' @param input_view A character name of the fst file that contains feature values
#' @param active_x feature name
#' @param readout A character of the score to be predicted ('bliss' or 'aoc')
#'
summarize_shap_values = function(feature_annotation, shap_view, active_x, readout,
                                 input_view) {

  # check if the queried data was cached before and return
  cache_vars = list(feature_annotation, shap_view, active_x, readout, input_view,
                    "shap")
  cache_hash = digest::digest(cache_vars)
  if (!dir.exists("data_caches"))
    dir.create("data_caches")

  cache_file = glue::glue("data_caches/{cache_hash}.RDS")
  if (file.exists(cache_file))
    return(readRDS(cache_file))

  # Get column names and select relevant ones
  columns = load_fst(fst_name = shap_view, peek = TRUE) %>%
    base::colnames()
  
  returns = list(df=NULL, x=NULL, y=NULL, color=NULL, fill=NULL)
  
  columns = c(base::grep(pattern = glue::glue("^{readout}"), x = columns, value=T),
              base::grep(pattern = "^\\.(meta|iden)", x = columns, value=T),
              base::grep(pattern = "^\\.response", x = columns, value=T),
              "feature")
  
  if (is.null(active_x)) {
    # summarize shap values by predefined top level feature categories
    active_x = "category"

    categories =  feature_annotation %>%
      dplyr::select(feature, !!active_x) %>%
      stats::na.omit()

    res = load_fst(fst_name = shap_view) %>%
      dplyr::filter(feature != "base_value") %>%
      dplyr::inner_join(categories, by="feature") %>%
      dplyr::group_by(.identifier_sample_name, !!rlang::sym(active_x)) %>%
      dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                          ~ abs(sum(.x, na.rm=TRUE))) %>%
      dplyr::group_by(!!rlang::sym(active_x)) %>%
      dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                          ~ mean(.x, na.rm=TRUE)) %>%
      tidyr::gather(-category, key="measure", value="shap") %>%
      tidyr::separate(measure, into=c("output", "fold"), sep="_y_") %>%
      dplyr::rename(x = !!active_x)

    returns$x_title = "All Feature Categories"

  } else {
    # check if the selected value is another category, otherwise check if a
    # feature was selected
    categories = tibble()
    # Lookup category
    if (all(active_x %in% unique(feature_annotation$category))) {
      categories =  feature_annotation %>%
        dplyr::filter(category %in% !!active_x) %>%
        dplyr::select(-gene, -category)
    } else if (all(active_x %in% unique(feature_annotation$molecular_biomarker))) {
      categories = feature_annotation %>%
        dplyr::filter(molecular_biomarker %in% active_x)
    } else if (all(active_x %in% unique(feature_annotation$chemical_structure))) {
      categories = feature_annotation %>%
        dplyr::filter(chemical_structure %in% active_x)
    }
    returns$x_title = stringr::str_to_title(gsub("_", " ", active_x))
    if (nrow(categories)  == 0) {
      
      res = load_fst(fst_name = shap_view,
                     filter   = list(feature = active_x),
                     cols     = columns)
      if (nrow(res) == 0) {
        return(NA)
      } else {
        
        inputs = load_fst(fst_name = input_view,
                          filter   = list(feature = active_x),
                          cols = c(".identifier_batch", ".identifier_sample_name", 
                                         ".metadata_moa_1", ".metadata_moa_2", 
                                         ".metadata_treatment_1", ".metadata_treatment_2", 
                                         "feature", "x"))
        if (is.numeric(inputs$x) & any(as.integer(inputs$x) != inputs$x)) {
          res = res %>%
            tidyr::gather(-c(dplyr::starts_with(c(".identifier", ".metadata",
                                                  ".response")),
                             feature, x),
                          key="measure", value="shap") %>%
            dplyr::filter(!is.na(shap)) %>%
            tidyr::separate(measure, into=c("output", "fold"), sep="_y_")
          
        } else {
          res = res %>%
            dplyr::left_join(inputs,  by = c(".identifier_batch", ".identifier_sample_name", 
                                             ".metadata_moa_1", ".metadata_moa_2", 
                                             ".metadata_treatment_1", ".metadata_treatment_2", 
                                             "feature")) %>%
            dplyr::filter(feature != "base_value") %>%
            dplyr::group_by(.identifier_sample_name, x) %>%
            dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                                ~ abs(sum(.x, na.rm=TRUE))) %>%
            dplyr::group_by(x) %>%
            dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                                ~ mean(.x, na.rm=TRUE)) %>%
            tidyr::gather(-x, key="measure", value="shap") %>%
            tidyr::separate(measure, into=c("output", "fold"), sep="_y_")
        }
      }
    } else {

      res = load_fst(fst_name = shap_view,
                     filter   = list(feature = categories$feature)) %>%
        dplyr::filter(feature != "base_value") %>%
        dplyr::inner_join(categories, by="feature")

      active_x = gsub(" ", "_", tolower(unique(active_x)))
      if (active_x %in% colnames(categories)) {

        res = res %>%
          dplyr::group_by(.identifier_sample_name, !!rlang::sym(active_x)) %>%
          dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                              ~ abs(sum(.x, na.rm=TRUE))) %>%
          dplyr::group_by(!!rlang::sym(active_x)) %>%
          dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                              ~ mean(.x, na.rm=TRUE)) %>%
          tidyr::gather(-!!active_x, key="measure", value="shap") %>%
          tidyr::separate(measure, into=c("output", "fold"), sep="_y_") %>%
          dplyr::rename(x = !!active_x)
      } else {
        res = res %>%
          dplyr::group_by(.identifier_sample_name, feature) %>%
          dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                              ~ abs(sum(.x, na.rm=TRUE))) %>%
          dplyr::group_by(feature) %>%
          dplyr::summarize_at(vars(dplyr::starts_with(readout)),
                              ~ mean(.x, na.rm=TRUE)) %>%
          tidyr::gather(-feature, key="measure", value="shap") %>%
          tidyr::separate(measure, into=c("output", "fold"), sep="_y_") %>%
          dplyr::rename(x = feature)
      }
    }
  }

  returns$df = res
  returns$x = "x"
  returns$y = "shap"
  returns$y_title = "Mean |SHAP values|"
  saveRDS(returns, cache_file)

  return(returns)
}