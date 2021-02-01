#' Encode drug names
#' @param x name of drug
#' @return numeric encode
#' @export

drug2idx <- function(x) {
    synergyddrapp::drug_summary %>%
        dplyr::select(drug_name) %>%
        dplyr::distinct() %>%
        dplyr::mutate(idx = dplyr::row_number()) %>%
        dplyr::filter(drug_name == !!x) %>%
        dplyr::pull(idx)
}

#' Encode mode of actions
#' @param x mode of action
#' @return numeric encode
#' @export
moa2idx <- function(x) {
    synergyddrapp::drug_summary %>%
        dplyr::select(`mode-of-action`) %>%
        dplyr::distinct() %>%
        dplyr::mutate(idx = dplyr::row_number()) %>%
        dplyr::filter(`mode-of-action` == !!x) %>%
        dplyr::pull(idx)
}
