#' Predict combination effect of a drug combination
#' @param 
#' @return
#' @export

 predict_combination_effect = function(X, feature_names, pred_target) {
    
     all_model_path = list.files(path = '../dependency/saved_models/', full.names=T,
                                 pattern = glue::glue("{pred_target}.*.model.txt$"))
     all_shap = c()

     preds = purrr::map(all_model_path, function(x) {
         reg = lightgbm::lgb.load(x)
         reg$predict(X)

         # TODO: shap score
         # lgb.trees <- lgb.model.dt.tree(reg) # First get a lgb tree
         # explainer <- lightgbmExplainer::buildExplainer(lgb.trees)
         # # compute contribution for each data point
         # pred.breakdown <- lightgbmExplainer::explainPredictions(reg, explainer, X)
     })
 }