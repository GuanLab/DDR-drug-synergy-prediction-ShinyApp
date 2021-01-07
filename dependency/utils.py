#!usr/bin/env/python3
#author: @rayezh
import pandas as pd
import numpy as np
from collections import defaultdict

class encode_categorical_feature():
    """ Categorical feature encoding object

    Parameters:
    -----------
    data: '../../feature/data/synergy_responses.tsv' pandas frame
    
    Yields:
    -------
    batch: dict
    cell_line: dict
    moa: dict
    remark: dict
    cancer_type: dict
    cancer_subtype: dict
    """

    data = pd.read_csv('./feature/QC/all_drugs_summary.csv')

    def __init__(self):
        self.moa = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(self.data['mode-of-action'].to_list())})
        self.treatment = defaultdict(lambda:np.nan, {j:i for i, j in enumerate(self.data['drug_name'].to_list())})

def boostrapping_confidence_interval(pred_all, ci):
    """ Boostrapping to get a 95 confidence interval for prediction performance
    
    Params
    ------
    pred_all: Numpy array
    ci: confidence interval
    
    Yields
    ------
    cor_mean: float 
        middle bound
    lb: float
        lower bound
    ub: float
        upper bound
    """
    import numpy as np
    import random
    # set random seed
    random.seed(0)
    
    
    # calculate overall correlation
    cor_mean = np.corrcoef(np.array(pred_all)[:,0], np.array(pred_all)[:,1])[0,1]
    print("Overall prediction/gold standard correlation is %.4f" % cor_mean)
    # start boostrapping ...
    cor_all = [] 
    for i in range(100):
        pred_new = random.choices(pred_all, k = len(pred_all))
        pred_new = np.array(pred_new)
        cor = np.corrcoef(pred_new[:,0], pred_new[:,1])[0,1]
        cor_all.append(cor)
    cor_all = sorted(cor_all)

    lb = cor_all[round(100*(0.5-ci*0.5))]
    ub = cor_all[round(100*(0.5+ci*0.5))]
    print("%d%s Confidence interval is: (%.4f, %.4f)" % (int(ci*100), '%', lb, ub))

    return cor_mean, lb, ub





