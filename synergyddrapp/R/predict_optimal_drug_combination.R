

#' Predict best drug combinations based on patient's molecular readouts
#' 
#' @param mol_df data.frame with molecualr readouts from patient
#' @return best drug combinations based on prediction
#' 
#' @export
predict_optimal_drug_combination = function(mol_df, data_path){

	# we can't used python
    #reticulate::source_python(system.file("python", "main_prediction.py", package="synergyddrapp")) 
    
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
    
    features = c('moa', 'drug_name', 'molecular', 'target_gene', 'chemical_structure', 'geneset')
    features = build_feature_dataset(data=df, mol=mol_df, features=features, data_path=data_path)
    feature_names = features[[1]] 
    Test_X = features[[2]] 

    # predict efficacy:
    auc = predict_combination_effect(X=Test_X, feature_names=feature_names, pred_target='aoc', data_path=data_path)
    pred_aoc = auc[[1]] 
    shap_aoc = auc[[2]]

    bliss = predict_combination_effect(Test_X, feature_names=feature_names, pred_target='bliss', data_path=data_path)
    pred_bliss = bliss[[1]] 
    shao_bliss = bliss[[2]]

    df$predicted_aoc = as.numeric(pred_aoc)
    df$predicted_bliss = as.numeric(pred_bliss)
    df
}
