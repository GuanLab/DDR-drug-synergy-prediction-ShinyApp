#!usr/bin/python3
#author:@rayezh
import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.ensemble import RandomForestRegressor

def train_RandomForest_model(X, Y):
    """ Train Random Forest models
    """
    regressor = RandomForestRegressor(n_estimators=500, random_state=0)
    regressor.fit(X, Y)
    return regressor

def train_LightGBM_model(X,Y, feature_name):
    """ Train LightGBM models
    
    Paramteres:
    -----------
    X: Numpy array
    Y: Numpy array
    feature_name: list of strings
    
    Yields:
    -------
    regressor: the light GBM model
    """
    param = {'boosting_type': 'gbdt',
            'objective': 'regression',
            'num_leaves': 20,
            'learning_rate': 0.05,
            'verbose': 0,
            'n_estimators': 1000,
            'reg_alpha': 2.0,
                   }
    # exclude categorical features that are not in feature set
    categorical_feature = ["Synergy_batch",
                "Cell_line",
                "Cancer_type",
                "Cancer_subtype",
                "Concentration",
                "Treatment_1_name",
                "Treatment_2_name",
                "Treatment_1_moa",
                "Treatment_2_moa"]
    
    categorical_feature = [x for x in categorical_feature if x in feature_name]

    # model training 
    train_data = lgb.Dataset(data = X,
            label = Y,
            feature_name = feature_name,
            categorical_feature = categorical_feature)
    regressor = lgb.train(param, train_data, num_boost_round=1000)
    return regressor
