#!usr/bin/python3
#author:@rayezh
import pandas as pd
import numpy as np
from tqdm import tqdm
import lightgbm as lgb
import shap
import matplotlib.pyplot as plt
from glob import glob
import pickle, json, os

from build_feature_dataset import *
from utils import *
from shap_analysis import *
from models import *

def train(train_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path):
    """ Train LightGBM synergy prediction model
    
    Params
    ------
    train_path: str
        path to training set
    target: str
        prediction target(aoc or bliss)
    model_path: str
        path where to save the trained  model

    Yields
    ------
    cor: float
        correlation between prediction and gold standard on training dataset
    """
    print("Loading Training dataset ...")
    Train = pd.read_csv(train_path, sep = '\t')
        
    # Exclude datasets with missing target score
    Train = Train[Train['.response_'+target].notna()]
        
    # Build Training feature set
    f_name, Train_X,Train_Y = build_feature_dataset(Train, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, if_train = True)
        
    print('start training model ...')
    regressor = train_LightGBM_model(Train_X, Train_Y, f_name)
    print('saving model ...')
    pickle.dump(regressor, open(model_path, 'wb'))
        
    beta_Y = regressor.predict(Train_X)
    cor = np.corrcoef(beta_Y, Train_Y)[0,1]
    print("Prediction-gold standard's correlation on training set: ", cor)
    
    return cor

def predict(test_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path):
    """ Perform trained models on test dataset
    
    Params
    ------
    test_path: str
        path to test set
    target: str
        prediction target(aoc or bliss)
    model_path: str
        path to saved model

    Yields
    ------
    cor: float
        correlation between prediction and gold standard on training dataset
    pred: Numpy array
        1st column: gold standard
        2nd column: prediction
    """
    print("Loading Test dataset ...")
    Test = pd.read_csv(test_path, sep = '\t')

    # Exclude datasets with missing target score
    Test = Test[Test['.response_'+target].notna()]
    if Test.shape[0] > 0:
        f_name, Test_X, Test_Y = build_feature_dataset(Test, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, if_train = False)
        print('Loading saved model ...')
        regressor = pickle.load(open(model_path, 'rb'))
        
        print('start prediction ...')
        pred_Y = regressor.predict(Test_X)
        
        # Save prediction and groud truth
        pred = np.array([Test_Y,pred_Y]).T
        
        cor = np.corrcoef(pred_Y, Test_Y)[0,1]
        print("Prediction-gold standard's correlation on test set:",cor)
    else:
        print("Empty Test Set!")
        cor = np.nan
        pred = np.nan
    return cor, pred
        
def specific_SHAP_analysis(data_path, subset_type, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features):
    """ carry out subset-specific SHAP analysis
    
    Params
    ------
    data_path: str
    subset_type: str
        the subset division used to carry out SHAP analysis. 
        for example: 
            'moa': mode-of-action specific
            'tissue': tissue specific
            'cell_line': cell line specific
    target: str
        predicted score ('aoc' or 'bliss')

    dep_gene: str
        gene to make dependency plot
    
    Yields
    ------
    """
    from glob import glob
    
    all_path = sorted(list(glob(data_path+'/fold_*')))

    for path in all_path:
        
        idx = path.split('_')[-1]
        
        # load trained model parameter
        model_path = 'monotherapy_'+target+'_'+str(idx)+'.model'
        
        # define path of subsets
        all_subset_path = glob(path+'/'+subset_type+'/*.tsv')
        assert len(all_subset_path)>0, 'File path: '+path+'/'+subset_type+'/ does not exist!'
        os.makedirs('./'+subset_type, exist_ok = True)
        
        for path in all_subset_path:
            subset_name = path.split('/')[-1].split('.')[0]
            subset_outpath = './'+subset_type+'/'+subset_name
            os.makedirs(subset_outpath, exist_ok = True)

            results = open(subset_outpath+'/cor_'+target+'.txt', 'a')
                
            # Load moa-specific dataset for SHAP analysis
            print("Loading dataset for mode-of-action",subset_name)
            Test = pd.read_csv(path, sep = '\t')
            Test = Test[Test['.response_'+target].notna()]
            if Test.shape[0] == 0:
                continue
            
            # prediction on test set
            cor, pred = predict(path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path)
            print('correlation =', cor)
            results.write('%.4f\n'%(cor))

            # SHAP analysis on the test set
            shap_df, shap_fig = SHAP_analysis(path,exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features,  model_path)

            # Save SHAP fig and dataframes
            shap_fig.savefig(subset_outpath+'/importance_scatter_'+target+'_'+str(idx)+'_lgb.pdf', format='pdf', dpi=1200, bbox_inches='tight')
            shap_df.to_csv(subset_outpath+'/SHAP_'+target+'_'+str(idx)+'.csv', index = False)
            
            results.close()

def five_fold_cross_validation(data_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features):
    """ Five fold cross validation of synergy prediction model

    Params
    ------
    data_path: str
    exclude_synergy_batch: boolean
        use synergy batch as feature or not
    exclude_cell_line: boolean
    exclude_cancer_type: boolean
    target: str
        predicted score (aoc or bliss)
    features: a list of strings
        features used to construct feature set
    
    Yields
    ------
    eva_df: Pandas dataframe
        evaluation results from k-fold cv
    eva_all: a str
        evalution reaulta from all k-fold; with 95 CI
    """
    from glob import glob
    
    eva_df = {"fold":[],"Pearson's r":[]} # predictions from all k-fold
    pred_all = [] # predictions from all k-fold
    
    all_path = sorted(list(glob(data_path+'/fold_*')))
    for path in all_path:

        idx = path.split('_')[-1]
        
        # start model training 
        train_path = path+'/Train.tsv'
        model_path = 'monotherapy_'+target+'_'+str(idx)+'.model'
        cor = train(train_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path)

        # prediction on test set
        test_path = path+'/Test.tsv'
        cor, pred = predict(test_path, exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features, model_path)
        np.savetxt(target+'_pred_'+str(idx)+'.txt', pred)
        
        # include fold prediction into all prediction
        pred_all.extend(pred.tolist())

        # write out test performance
        eva_df["fold"].append(idx)
        eva_df["Pearson's r"].append(cor)
        
        # SHAP analysis on the overall test set
        shap_df, shap_fig = SHAP_analysis(test_path,exclude_synergy_batch, exclude_cell_line, exclude_cancer_type, target, features,  model_path)

        # Save SHAP fig and dataframes
        shap_fig.savefig('importance_scatter_'+target+'_'+str(idx)+'_lgb.pdf', format='pdf', dpi=1200, bbox_inches='tight')
        shap_df.to_csv('SHAP_'+target+'_'+str(idx)+'.csv', index = False)

    eva_df = pd.DataFrame.from_dict(eva_df)

    # calculate overall prediction performance confidence interval
    pred_all =  np.array(pred_all)
    mb, lb, ub = boostrapping_confidence_interval(pred_all, 0.95)
    eva_all = "mean[95CI]: %.4f[%.4f, %.4f]" % (mb, lb, ub)
    
    return eva_df, eva_all

