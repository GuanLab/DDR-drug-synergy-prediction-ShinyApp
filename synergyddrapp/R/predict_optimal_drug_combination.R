

#' Predict best drug combinations based on patient's molecular readouts
#' 
#' @param mol_df data.frame with molecualr readouts from patient
#' @return best drug combinations based on prediction
#' 
#' @export
predict_optimal_drug_combination = function(mol_df, data_path) {

    reticulate::source_python(system.file("python", "main_prediction.py", package="synergyddrapp")) 
    
    all_drugs = synergyddrapp::drug_summary

    # create combination of all drugs and MoAs
    .combn_pair = function(x, col) {
        combn(x[[col]], m=2, simplify=F) %>% 
            unlist %>% 
            matrix(., ncol=2, byrow=T) %>%
            data.frame()
    }
    tratment_comb_2 = .combn_pair(all_drugs, "drug_name") %>%
        dplyr::rename(.metadata_treatment_1=X1, .metadata_treatment_2=X2)
    moa_comb_2 = .combn_pair(all_drugs, "mode-of-action") %>%
        dplyr::rename(.metadata_moa_1=X1, .metadata_moa_2=X2)
    
    df = cbind(moa_comb_2, tratment_comb_2)
    
    features = synergyddrapp::build_feature_dataset(data=df, mol=mol_df,  target="aoc", data_path=data_path)

    # predict efficacy:
    pred_aoc = synergyddrapp::predict_combination_effect(X=features, pred_target='aoc')

    pred_bliss = synergyddrapp::predict_combination_effect(X=features, pred_target='bliss')

    df %>%
        dplyr::mutate(predicted_aoc = pred_aoc,
                      predicted_bliss = pred_bliss)
}
