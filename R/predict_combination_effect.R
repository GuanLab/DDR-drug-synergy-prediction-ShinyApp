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

  predict_with_fallback = function(reg, newdata, predcontrib) {
    # Some serialized model wrappers call an internal predictor with
    # a stale reshape argument. Fall back to stable prediction APIs.
    primary = tryCatch({
      reg$predict(newdata, predcontrib = predcontrib)
    }, error = function(e) e)

    if (!inherits(primary, "error"))
      return(primary)

    use_fallback = grepl("unused argument \\(reshape = reshape\\)",
                         primary$message, fixed = FALSE)
    if (!use_fallback)
      stop(primary)

    fallback_predict = function(model) {
      tryCatch({
        stats::predict(model, newdata = newdata, predcontrib = predcontrib)
      }, error = function(e) e)
    }

    candidates = list(reg)
    if (!is.null(reg$finalModel))
      candidates = c(candidates, list(reg$finalModel))
    if (!is.null(reg$booster))
      candidates = c(candidates, list(reg$booster))

    for (candidate in candidates) {
      res = fallback_predict(candidate)
      if (!inherits(res, "error"))
        return(res)
    }

    stop(primary)
  }

  preds = purrr::map(all_model_path, function(x) {
    reg = readRDS(x)
    predict_with_fallback(reg = reg,
                          newdata = as.matrix(X),
                          predcontrib = predcontrib)
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
