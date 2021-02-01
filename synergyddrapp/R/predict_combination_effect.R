#' @export
predict_combination_effect = function(X, pred_target) {
    
    all_model_path = list.files(path = system.file("extdata", "saved_models", package="synergyddrapp"), 
                                full.names=T, pattern = glue::glue("{pred_target}.*.model.txt$"))
    all_shap = c()
    
    preds = purrr::map(all_model_path, function(x) {
        reg = lightgbm::lgb.load(x)
        reg$predict(as.matrix(X))
    })

    data.frame(t(matrix(unlist(preds), nrow=length(preds), byrow=T))) %>%
        rowMeans(.)
}