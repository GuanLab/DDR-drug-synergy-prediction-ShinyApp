
#' Summary function for input values
#'
#' This function summarizes the input values of the test data for a defined
#' set of features
#'
#' @param feature_annotation A data.frame defining categories of features 
#' (see [create_feature_categories()])
#' @param view A character name of the view that contains input values
#' @param active_x A character feature name for that the input statistics should
#' be summarized
#'
summarize_input_values = function(feature_annotation, view, active_x) {

  # check if the queried data was cached before and return
  cache_vars = list(feature_annotation, view, active_x, "input")
  cache_hash = digest::digest(cache_vars)
  if (!dir.exists("data_caches"))
    dir.create("data_caches")

  cache_file = glue::glue("data_caches/{cache_hash}.RDS")
  if (file.exists(cache_file))
    return(readRDS(cache_file))
  
  # Get column names and select relevant ones
  columns = load_fst(fst_name = view, peek = TRUE) %>%
    base::colnames()
  
  returns = list(df=NULL, x=NULL, y=NULL, color=NULL, fill=NULL)
  
  columns = c(base::grep(pattern = "^\\.(meta|iden)", x = columns, value=T),
              base::grep(pattern = "^\\.response", x = columns, value=T),
              "feature", "x")
  
  if (is.null(active_x)) {
    # summarize shap values by predifined top level feature categories

    returns$df =  feature_annotation %>%
      dplyr::select(feature, category, molecular_biomarker,
                    chemical_structure) %>%
      dplyr::distinct() %>%
      tidyr::unite(col = "feature_class", -feature, na.rm=TRUE) %>%
      dplyr::group_by(feature_class) %>%
      dplyr::summarise(x = n()) %>%
      tidyr::separate(feature_class, into= c("category", "subcategory"),
                      sep="_") %>%
      dplyr::mutate(subcategory = ifelse(is.na(subcategory),
                                         category, subcategory)) %>%
      dplyr::group_by(category) %>%
      dplyr::mutate(Total=sum(x)) %>%
      dplyr::arrange(dplyr::desc(Total), dplyr::desc(x)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(category = forcats::fct_inorder(factor(category)),
                    subcategory = forcats::fct_inorder(factor(subcategory)))

    returns$x    = "category"
    returns$y    = "x"
    returns$fill = "subcategory"
    returns$y_title = "Number of Features"
    returns$x_title = "All Feature Categories"
    returns$legend_title = "Subcategory"

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

    if (nrow(categories) > 0) {
      categories = categories %>%
        dplyr::select_if(~any(!is.na(.x)))

      res = load_fst(fst_name = view, cols = columns,
                     filter   = list(feature = categories$feature)) %>%
        dplyr::inner_join(categories, by="feature")

      active_x = gsub(" ", "_", tolower(unique(active_x)))
      if (active_x %in% colnames(categories)) {

        returns$df = res %>%
          dplyr::distinct(!!rlang::sym(active_x), feature) %>%
          dplyr::count(!!rlang::sym(active_x)) %>%
          dplyr::arrange(dplyr::desc(n)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(!!active_x := forcats::fct_inorder(factor(!!rlang::sym(active_x))))

        returns$x = active_x
        returns$x_title = stringr::str_to_title(gsub("_", " ", active_x))
        returns$fill = NULL
        returns$color = NULL
        returns$y = "n"
        returns$y_title = "Number of Features"


      } else {

        if (is.integer(res$x) | all(res$x == as.integer(res$x))) {

          feature_variance = res %>%
            dplyr::group_by(feature) %>%
            dplyr::summarise(median_value = median(x)) %>%
            dplyr::arrange(-median_value)

          returns$df = res %>%
            dplyr::group_by(feature, x) %>%
            dplyr::summarise(count = n()) %>%
            dplyr::inner_join(feature_variance) %>%
            dplyr::mutate(feature = factor(feature, levels=feature_variance$feature))

          if (!class(returns$df[["x"]]) %in% c("character", "factor"))
            returns$df = returns$df %>%
             dplyr::mutate(x = factor(x, levels=as.character(unique(.$x)[order(unique(.$x))])))

          returns$x = "feature"
          returns$y = "count"
          returns$y_title = "Number of Features"
          returns$x_title = stringr::str_to_title(gsub("_", " ", active_x))
          returns$fill = "x"
          returns$legend_title = "Value"
        } else {

          feature_variance = res %>%
            dplyr::group_by(feature) %>%
            dplyr::summarise(var = median(x)) %>%
            dplyr::arrange(-var)

          returns$df = res %>%
            dplyr::inner_join(feature_variance) %>%
            dplyr::group_by(feature) %>%
            dplyr::summarise(min=min(x, na.rm = TRUE),
                             max=max(x, na.rm = TRUE),
                             middle=median(x, na.rm = TRUE),
                             lower = quantile(x, .25, na.rm = TRUE),
                             upper = quantile(x, .75, na.rm = TRUE))%>%
            dplyr::mutate(feature = factor(feature, levels=feature_variance$feature))
          returns$x = "feature"
          returns$x_title = stringr::str_to_title(gsub("_", " ", active_x))
        }
      }
    }

  }
  saveRDS(returns, cache_file)

  return(returns)
}
