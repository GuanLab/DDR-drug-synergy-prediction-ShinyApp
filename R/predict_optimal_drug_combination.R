#' Predict best drug combinations based on patient's molecular readouts
#' 
#' @param mol_df data.frame with molecualr readouts from patient
#' @return best drug combinations based on prediction
#' 
#' @export
predict_optimal_drug_combination = function(mol_df) {

    
    tratment_comb_2 <- .combn_pair(synergyddrapp::drug_summary,
                                   "drug_name") %>%
        dplyr::rename(.metadata_treatment_1=X1, .metadata_treatment_2=X2)
    
    moa_comb_2 <- .combn_pair(synergyddrapp::drug_summary, "mode-of-action") %>%
        dplyr::rename(.metadata_moa_1=X1, .metadata_moa_2=X2)
    
    df <- cbind(moa_comb_2, tratment_comb_2)
    
    # create a set of features
    features <- build_feature_dataset(data=df, mol=mol_df, target="aoc")

    # predict efficacy:
    pred_aoc         <- predict_combination_effect(X=features, pred_target='aoc')
    feature_shap_aoc <- predict_combination_effect(X=features, pred_target='aoc',
                                                  predcontrib=TRUE)
     
    # create a set of features
    features <- build_feature_dataset(data=df, mol=mol_df, target="bliss")
    
    pred_bliss         <- predict_combination_effect(X=features, 
                                                    pred_target='bliss')
    feature_shap_bliss <- predict_combination_effect(X=features, 
                                                    pred_target='bliss', 
                                                    predcontrib=TRUE)

    # append predictions to treatment/moa combinations
    prediction <- df %>%
        dplyr::mutate(predicted_aoc = pred_aoc,
                      predicted_bliss = pred_bliss)

    list(predictions=prediction, bliss_contributions=feature_shap_bliss, 
         aoc_contribution=feature_shap_aoc)
}
