#' @export
predict_combination_effect = function(X, pred_target, predcontrib=FALSE) {
    
    all_model_path = list.files(path = system.file("extdata", "saved_models", package="synergyddrapp"), 
                                full.names=T, pattern = glue::glue("{pred_target}.*.model.txt$"))
    all_shap = c()
    
    preds = purrr::map(all_model_path, function(x) {
        reg = lightgbm::lgb.load(x)
        reg$predict(as.matrix(X),  predcontrib=predcontrib)
    })

    if (!predcontrib) {
        data.frame(t(matrix(unlist(preds), nrow=length(preds), byrow=T))) %>%
            rowMeans(.)
    } else {
       purrr::map(preds, function(x) {
           x = dplyr::as_tibble(x[, -ncol(x)])
           colnames(x) = colnames(X)
           x
           }) 
           
    }
}