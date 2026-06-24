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

  data = as.matrix(X)

  preds = purrr::map(all_model_path, function(x) {
    reg = readRDS(x)
    predictor = reg$to_predictor()
    predictor_formals = tryCatch(names(formals(predictor$predict)),
                                 error = function(...) character())
    predictor_args = list(
      data = data,
      start_iteration = 0L,
      num_iteration = reg$best_iter,
      rawscore = FALSE,
      predleaf = FALSE,
      predcontrib = predcontrib,
      header = FALSE
    )

    if ("reshape" %in% predictor_formals)
      predictor_args$reshape = FALSE

    do.call(predictor$predict, predictor_args)
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
