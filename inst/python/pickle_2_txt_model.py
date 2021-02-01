import pickle
from glob import glob
for pred_target in ['aoc', 'bliss']:
    all_model_path = glob('../saved_models/*_'+pred_target+'_*.model')
    for model_path in all_model_path:
        print("loading "+pred_target+" model...")
        regressor = pickle.load(open(model_path, 'rb'))
        regressor.save_model(model_path+'.txt', num_iteration=regressor.best_iteration)
