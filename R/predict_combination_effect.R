#' Predict treatment combination effect
#'
#' @param X A matrix of all features as used for training
#' @param pred_target A character for the score to be predicted ('aoc' or bliss)
#' @param predcontrib A logical for returning pe-feature contributions or not
#' @param path A character of the path to model files
#' @param threads The number of threads to be used by \link{lightgbm}
#' (name: '{pred_target}_{cnt}.model.txt')
#'
predict_combination_effect = function(X, pred_target, predcontrib=FALSE,
                                      path=NULL, threads=1) {


  all_model_path = list.files(path = path,
                              full.names=T,
                              pattern = pred_target)

  preds = purrr::map(all_model_path, function(x) {
    reg = lightgbm::readRDS.lgb.Booster(x)
    reg$predict(as.matrix(X),  predcontrib=predcontrib, num_threads=threads)
  })

  all_shap = c()
  
  if (!predcontrib) {
    data.frame(t(matrix(unlist(preds), nrow=length(preds), byrow=T))) %>%
      rowMeans()
  } else {
    purrr::map(preds, function(x) {
      colnames(x) = c(colnames(X), "base_value")
      dplyr::as_tibble(x)
    })
  }
}
