import pandas as pd
import numpy as np
import lightgbm as lgb
import shap
from glob import glob
import json, os, re, pickle
from collections import defaultdict
from tqdm import tqdm

class encode_categorical_feature():
    """ Categorical feature encoding object

    Parameters:
    -----------
    data: 
    
    Yields:
    -------
    batch: dict
    cell_line: dict
    moa: dict
    remark: dict
    cancer_type: dict
    cancer_subtype: dict
    """

    data = pd.read_csv('./dependency/feature/QC/all_drugs_summary.csv')

    def __init__(self):
        self.moa = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(self.data['mode-of-action'].to_list())})
        self.treatment = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(self.data['drug_name'].to_list())})

def build_feature_dataset(data, mol, features):
    """ Construc feature set
    
    Parameters:
    -----------
    
    data: Pandas dataframe
        synergy response data in training or testing dataset
    
    exclude_synergy_batch: boolean
        if True, exclude synergy batch.
        Used for cross-batch validation

    exclude_cell_line: boolean
        if True, exclude cell line
        Used for cross-cell line validation

    exclude_cancer: boolean
        if True, exclude cancer tissue
        Used for cross-indication validation

    target: string
        prediction target('aoc' or 'bliss')
    
    features: list
        features to use in prediction.
        monotherapy response
        molecluar
        chemical_structure
        crispr

    if_train: boolean
        if construct training data, true; if test data, false
    
    return_label: boolean
        option to return label (for prediction)
        or not to return label (for held-out set prediction, where the label is unknown))
    
    Yields:
    -------
    feature_name: list of strings
        feature names
    X: Numpy array
        feature matrix
    Y: Numpy array
        predicted label(aoc score or bliss score)
    """

    # load feature datasets
    all_chemical_structure = pd.read_csv('./dependency/feature/chemical_structure_features/all_chemical_structure.csv', index_col = 0)
    geneset = pd.read_csv('./dependency/feature/geneset_features/geneset_features.csv', index_col = 0)
    ceres_cancer_dependency = json.load(open('./dependency/feature/cancer_dependency_features/integrated_Sanger_Broad_essentiality_matrices_20200402/CERES_FC_dep.json', 'r'))
    crisprclear_cancer_dependency = json.load(open('./dependency/feature/cancer_dependency_features/integrated_Sanger_Broad_essentiality_matrices_20200402/CRISPRcleanR_FC_dep.json', 'r'))

    # load drug-specific target gene information
    drug2gene = json.load(open('./dependency/feature/target_gene/drug2gene.json','r'))
    drug2gene = defaultdict(lambda:[], drug2gene)

    # load gene network information from mousenet
    network = pickle.load(open('./dependency/feature/mousenet_features/target_gene_network.pkl','rb'))

    # categotical encoding
    encode = encode_categorical_feature() 

    def build_features(row, mol, features):
        """ Build up all types of features
        
        Params
        ------
        row: Pandas row series
            a single line of synergy experimental data
        features: a list of str
            features used in feature set
        exclude_synergy_batch: boolean
        exclude_cell_line: boolean
        exclude_cancer_type: boolean

        Yields:
        -------
        feature_names: list
        feature_values: list
        """

        def moa_feature(row):
            """ Mode-of-action as categorical feature

            Params
            ------
            row: Pandas row series

            Yields
            ------
            feature_names: list
            feature_values: list
            """
            feature_names = []
            feature_values = []
            for i in ['1','2']:
                moa = row['.metadata_moa_'+i]
                feature_names.extend(['Treatment_'+i+'_moa'])
                feature_values.extend([encode.moa[moa]])

            return feature_names, feature_values

        def drug_name_feature(row):
            """ Drug name as categorical feature
            
            Params
            ------
            row:Pandas row series

            Yields
            ------
            feature_names: list
            feature_values: list
            """
            feature_names = []
            feature_values = []
            for i in ['1','2']:
                treatment = row['.metadata_treatment_'+i]
                feature_names.extend(['Treatment_'+i+'_name'])
                feature_values.extend([encode.treatment[treatment]])

            return feature_names, feature_values

        def molecular_feature(row, mol, include_target_gene = False, include_network =  False):
            """ Molecular features
            
            Params
            ------
            row: Pandas row series
            include_target_gene: boolean
                if include target gene information or not
            include_network: boolean
                include network information or not
                can only be set to true if include_target_gene is true

            Yields
            ------
            feature_names: list
            feature_values: list
            """

            mol = mol.iloc[0].copy()
            
            # include target gene information
            if include_target_gene:
                # get target genes from treatment 1 and 2
                target_genes = set(drug2gene[row['.metadata_treatment_1']]+drug2gene[row['.metadata_treatment_2']])
                target_genes_exp = [g+'_exp' for g in target_genes]
                # set <target gene>_exp to 0
                for exp in target_genes_exp:
                    if exp in mol.index:
                        mol[exp] = 0
                    else:
                        pass

                if include_network:
                    # maximum relationship with target genes
                    p2target = defaultdict(lambda:0)
                    for target_gene in target_genes:
                        for k, v in network[target_gene].items():
                            if v>p2target[k]:
                                p2target[k] = v
                    for k,v in p2target.items():
                        if k+'_exp' in mol.index:
                            mol[k+'_exp'] = mol[k+'_exp']*(1-v)

            feature_names = list(mol.index)
            feature_values = list(mol)
            
            return feature_names, feature_values

        def chemical_structure_feature(row):
            """ Drug features:
                1. treatment name (categorical feature)
                2. drug chemical structure feature
            
            Params
            ------
            row: Pandas row series

            Yields
            ------
            feature_names: list
            feature_values: list
            """
            feature_names = []
            feature_values = []
            for i in ['1','2']:
                treatment = row['.metadata_treatment_'+i]
                # extend feature name
                feature_names.extend(['Treatment_'+i+'_'+c for c in all_chemical_structure.columns.to_list()])
                # extend feature values
                if treatment in all_chemical_structure.index:
                    feature_values.extend(all_chemical_structure.loc[treatment].to_list())
                else:
                    feature_values.extend([np.nan]*len(all_chemical_structure.columns))
                    if treatment == "GAMMA":
                        pass
                    else:
                        print("No chemical structure data of "+treatment+"! Replace with NaN!")
            
            return feature_names, feature_values
        
        def geneset_feature(row):
            """ Geneset feature 
            
            Params
            ------

            Yields
            ------
            feature_value: a vector of n genesets by 1.
                for each geneset, it returns a number denoting how many genes targeted by treatments are in the pathway
            """
            feature_names = list(geneset.index)
            feature_values = np.zeros(len(feature_names))

            # get target genes from treatment 1 and 2
            target_genes = set(drug2gene[row['.metadata_treatment_1']]+drug2gene[row['.metadata_treatment_2']])

            feature_names = ["Geneset_"+x for x in list(geneset.index)]
            for g in target_genes:
                try:
                    feature_values += np.array(geneset[g])
                except:
                    pass
                    #print("no geneset feature!")
            feature_values = list(feature_values)
            
            return feature_names, feature_values

        def cancer_dependency_feature(row):
            """ Cancer dependency feature
            
            Params
            ------
            row: Pandas row series
            
            Yields
            ------
            feature_names: list
            feature_values: list
            """
            cell = row['.identifier_sample_name']
            moas = [row['.metadata_moa_1'], row['.metadata_moa_2']]
            
            feature_names = []
            feature_values = []

            for i in ['ceres', 'crispr']:
                f_name = 'CRISPR_cancer_dependency_'+i
                if i == 'ceres':
                    df = ceres_cancer_dependency
                if i == 'crispr':
                    df = crisprclear_cancer_dependency
                for j in ['max']: # max dependency
                    f_name = f_name+'_'+j
                    feature_names.append(f_name)
                    try:
                        v = []
                        for moa in moas:
                            v.append(df[moa][cell][j])
                        if j == 'max':
                            v = max(v)
                        elif j == 'mean':
                            v = np.mean(v)
                        elif j == 'min':
                            v = min(v)
                    except:
                        v = np.nan
                    feature_values.append(v)

            return feature_names, feature_values

        # aggregate all features above:
        feature_names = []
        feature_values = []
        
        if 'moa' in features:
            moa_name, moa_feature = moa_feature(row)
            feature_names += moa_name
            feature_values += moa_feature
        
        if 'drug_name' in features:
            dr_name, dr_feature = drug_name_feature(row)
            feature_names += dr_name
            feature_values += dr_feature

        if 'monotherapy' in features:
            mono_name, mono_feature = monotherapy_feature(row)
            feature_names += mono_name
            feature_values += mono_feature

        if 'molecular' in features:
            include_target_gene = False
            include_network = False
            if 'target_gene' in features:
                include_target_gene = True
                if 'network' in features:
                     include_network = True
            mol_name, mol_feature = molecular_feature(row, mol, include_target_gene, include_network)
            
            feature_names += mol_name
            feature_values += mol_feature

        if 'chemical_structure' in features:
            ch_name, ch_feature = chemical_structure_feature(row)
            feature_names += ch_name
            feature_values += ch_feature
        
        if 'geneset' in features:
            g_name, g_feature = geneset_feature(row)
            feature_names += g_name
            feature_values += g_feature

        if 'dependency' in features:
            dep_name, dep_feature = cancer_dependency_feature(row)
            feature_names += dep_name
            feature_values += dep_feature
        
        return feature_names, feature_values
    
    # initiat feature set(X) and label set(Y)
    X = []
    feature_names = []
    # Read by row, construct total feature sets from synergetic experiment information
    data = data.reset_index(drop = True)
    for idx, row in tqdm(data.iterrows(), total = data.shape[0]):
        # feature name list
        feature_names, feature_values = build_features(row, mol, features)
        if idx == 0:
            print("Total features:",len(feature_names))

        X.append(feature_values)
        
    X = np.array(X)
    return feature_names, X

def predict_combination_effect(X,feature_names,pred_target):
    """ Predict effect of drug combinations
    
    Params
    ------
    pred_target: str
        aoc or bliss
    Yields
    ------
    """
    all_model_path = glob('./dependency/saved_models/*_'+pred_target+'_*.model')
    all_pred = []
    all_shap = []
    for model_path in all_model_path:
        print("loading "+pred_target+" model...")
        regressor = pickle.load(open(model_path, 'rb'))

        # make prediction
        pred_Y = regressor.predict(X)
        all_pred.append(pred_Y)

        # shap analysis
        shap_values = shap.TreeExplainer(regressor).shap_values(X)
        all_shap.append(shap_values)
    all_pred = np.mean(np.array(all_pred),axis=0)
    shap_df = pd.DataFrame(np.array(shap_values), columns = feature_names)
    return all_pred, shap_df

def predict_optimal_drug_combination(mol_df):
    """ Predict best drug combinations based on patient's molecular readouts
    
    Params
    ------
    mol_df: Pandas dataframe
        molecualr readouts from patient
    
    Yields
    ------
    best drug combinations based on prediction

    """
    all_drugs =pd.read_csv('./dependency/feature/QC/all_drugs_summary.csv')

    # create a list all drug combinations as test set
    df = {'.metadata_moa_1':[],'.metadata_moa_2':[],'.metadata_treatment_1':[],'.metadata_treatment_2':[]}
    for i in range(all_drugs.shape[0]-1):
        for j in range(i+1, all_drugs.shape[0]):
            #if all_drugs['mode-of-action'].iloc[i] in ['ATMi' 'ATRi', 'DNAPKi']:
            df['.metadata_treatment_1'].append(all_drugs['drug_name'].iloc[i])
            df['.metadata_moa_1'].append(all_drugs['mode-of-action'].iloc[i])
            df['.metadata_treatment_2'].append(all_drugs['drug_name'].iloc[j])
            df['.metadata_moa_2'].append(all_drugs['mode-of-action'].iloc[j])

    df = pd.DataFrame.from_dict(df)
    features = ['moa', 'drug_name', 'molecular', 'target_gene', 'chemical_structure', 'geneset']
    feature_names, Test_X = build_feature_dataset(df, mol_df, features)

    # predict efficacy:
    pred_aoc, shap_aoc = predict_combination_effect(Test_X, feature_names, 'aoc')
    pred_bliss, shao_bliss = predict_combination_effect(Test_X,feature_names, 'bliss')
    df['predicted_aoc'] = list(pred_aoc)
    df['predicted_bliss'] = list(pred_bliss)
    return df
