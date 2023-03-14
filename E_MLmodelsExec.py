# import necessary libraries
from random import seed
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import f1_score
from sklearn.metrics import r2_score
from sklearn import model_selection
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, recall_score, confusion_matrix
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import glob, os
from tqdm import tqdm
from numpy import mean
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.utils import resample
import random
import numpy as np
from sklearn import datasets, metrics, model_selection, svm
from lightgbm import LGBMClassifier

import lightgbm as lgb
from sklearn.model_selection import cross_val_score
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
from hyperopt.pyll import scope
import time
import warnings
warnings.filterwarnings('ignore')

# Importing core libraries
import numpy as np
import pandas as pd
from time import time
import pprint
import joblib
from functools import partial
from bayes_opt import BayesianOptimization
from sklearn.model_selection import StratifiedKFold as skf

# Suppressing warnings because of skopt verbosity
import warnings
warnings.filterwarnings("ignore")

# Classifiers
import lightgbm as lgb

# Model selection
from sklearn.model_selection import KFold, StratifiedKFold

# Metrics
from sklearn.metrics import mean_squared_error
from sklearn.metrics import make_scorer

# Skopt functions
from skopt import BayesSearchCV
from skopt.callbacks import DeadlineStopper, DeltaYStopper
from skopt.space import Real, Categorical, Integer

# Data processing
from sklearn.preprocessing import OrdinalEncoder
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer

# set seed for reproducibility
os.environ["PL_GLOBAL_SEED"] = str(seed)
random.seed(seed)

# define function to get list of models to evaluate
def get_models():
	models = list()
	models.append(LGBMClassifier(max_depth=10, n_estimators=500))
	return models

# Get the current working directory
cwd = os.getcwd()

# Construct paths to input and output files relative to the current working directory
input_file_path = os.path.join(cwd, 'Python\\libFiles\\downsample')
output_file_path = os.path.join(cwd, 'Python\\libFiles\\processed')
input_kmer_path = os.path.join(cwd, 'Dataset')
output_results_path = os.path.join(cwd, 'Resultados')

# load list of selected samples
antibiotic_dfs ={}
antibiotics = ['ceftazidime', 'ciprofloxacin','meropenem', 'tobramycin']
for antibiotic in antibiotics:
    antibiotic_dfs[antibiotic] = pd.read_csv(os.path.join(input_file_path,antibiotic+'_AMR.csv'), index_col=False, sep=';', header = 0)
    #antibiotic_dfs[antibiotic] = pd.read_csv(os.path.join(cwd,antibiotic+'_AMR.csv'), index_col=False, sep=';', header = 0)
    print(antibiotic, antibiotic_dfs[antibiotic].shape)
    #antibiotic_dfs[antibiotic].set_index('Index', inplace=True, drop=True)

# call function to get list of models
models = get_models()

# initialize KFold cross-validation
cont=0
cv = KFold(n_splits=10, shuffle=True, random_state=1)



############################# inicio ###################################



############################# Escolha do K ###################################

# Iterates over the values of k to read the corresponding k-mer feature dataset
for i in tqdm(range(6,7)):
    # Reads the k-mer feature dataset file
    

    # Iterates over the models list
    for model in models:
        # Prints the name of the current model being used
        print('\n\n'+type(model).__name__)
        
        # Opens the metadata file for the current model and delta values to store the results
        with open(os.path.join(output_results_path,type(model).__name__ +'.txt'), 'a+') as f:
        #with open(os.path.join(cwd,type(model).__name__ +'.txt'), 'a+') as f:
        #with open('/home/gabrielteixeirasousa/'+type(model).__name__ +str(delta)+'.txt', 'a+') as f:
        
            try: 
                # Checks if the metadata file already exists
                with open(os.path.join(output_results_path,type(model).__name__ +'.txt'), 'r') as f2:
                #with open(os.path.join(cwd,type(model).__name__ +'.txt'), 'r') as f2:
                #with open('/home/gabrielteixeirasousa/'+type(model).__name__ +str(delta)+'.txt', 'r') as f2:
                    # Reads the first line of the metadata file
                    primeiralinha = f2.readline()
                print(primeiralinha)
                # Checks if the first line of the metadata file is the expected header, if not, writes it
                if str(primeiralinha) !='index;ceftazidime;ciprofloxacin;meropenem;tobramycin;\n':
                    f.write('index;ceftazidime;ciprofloxacin;meropenem;tobramycin;')
            except: 
                # Writes the header in case the metadata file does not exist
                f.write('index;ceftazidime;ciprofloxacin;meropenem;tobramycin;')
            

            f.write('\n'+str(i)+';')


            for antibiotic in tqdm(antibiotic_dfs):
                #atual_XIS = filter_by_genome(xis, antibiotic)
                #xis = pd.read_csv(os.path.join(cwd,'kmer' +str(i) +'.csv'), sep=';', header = 0)
                xis = pd.read_csv(os.path.join(input_kmer_path,'kmer' +str(i) +'.csv'), sep=';', header = 0)
                #xis = pd.read_csv('/home/gabrielteixeirasousa/kmer' +str(i) +'.csv', sep=';', header = 0)
                # Sets the index of the DataFrame to the first column and drops it
                xis.set_index('Unnamed: 0', inplace=True, drop=True)
                xis.fillna(0, inplace=True)    


                for k in xis.iterrows():
                    arg1 = str(k[0]).split('.')[0]
                    if arg1 not in antibiotic_dfs[antibiotic]['Genome Name'].tolist():
                        xis = xis.drop(k[0])


                antibiotic_dfs[antibiotic].fillna(0, inplace=True)
                yps = antibiotic_dfs[antibiotic][antibiotic]
                X_train, X_test, y_train, y_test = model_selection.train_test_split(xis, yps ,test_size=0.2, 
                                                                                    random_state=0, stratify=antibiotic_dfs[antibiotic][antibiotic])

                '''print(antibiotic, y_train.value_counts())
                # Separate majority and minority classes
                df_ones = y_train[y_train==1]
                df_zeros = y_train[y_train==0]

                if len(df_ones)>= len(df_zeros):
                    df_majority = df_ones
                    df_minority = df_zeros
                    maj = 'one'
                else:
                    df_majority = df_zeros
                    df_minority = df_ones
                    maj = 'zero'
                
                # Calculate the proportion of the majority class
                proportion_majority = len(df_majority) / len(y_train)
                
                if maj == 'one' and proportion_majority < 0.6:
                    # Downsample majority class
                    #srt = math.sqrt(    4*(     (len(df_majority)/0.55)  - len(df)**2  )   )
                    
                    #x = int(( (-2*len(df) + 1/0.55) + srt )/2)
                    df_minority_downsampled = resample(df_minority, replace=False, n_samples=int(len(df_minority)*0.5), random_state=42)
                    y_train_df_downsampled = pd.concat([df_minority_downsampled, df_majority])
                    print(y_train_df_downsampled.value_counts())
                elif maj == 'zero' :
                    df_majority_downsampled = resample(df_majority, replace=False, n_samples=int(len(df_minority)*0.5), random_state=42)
                    y_train_df_downsampled = pd.concat([df_majority_downsampled, df_minority])

                    #df_downsampled = df
                    print(y_train_df_downsampled.value_counts())
                else:
                    y_train_df_downsampled = y_train


                y_train_df_downsampled = y_train_df_downsampled.iloc[y_train.index.values.tolist()]
                #y_train_df_downsampled = y_train_df_downsampled.sort_index(ascending=True)
                antibiotic_dfs[antibiotic] = antibiotic_dfs[antibiotic].loc[antibiotic_dfs[antibiotic].index.isin(y_train_df_downsampled.index)]

                # create a set of the genome names in df2
                genome_names = set(antibiotic_dfs[antibiotic]['Genome Name'])

                # filter df1 based on whether its index is in the set of genome names
                X_train_df_downsampled = X_train
                for k in X_train_df_downsampled.iterrows():
                    if k[0].split('.fna')[0] not in genome_names:
                        X_train_df_downsampled = X_train_df_downsampled.drop(k[0])

                X_train_df_downsampled = X_train_df_downsampled.sort_index(ascending=True)


                #X_train = X_train[X_train.index.values.tolist().isin(genome_names)]'''

                # define objective function
                def lgbm_cv(n_estimators, learning_rate, max_depth, num_leaves, min_child_samples):
                    params = {'n_estimators': int(round(n_estimators)),
                            'learning_rate': learning_rate,
                            'max_depth': int(round(max_depth)),
                            'num_leaves': int(round(num_leaves)),
                            'min_child_samples': int(round(min_child_samples))}
            
                    # define cross-validation splitter
                    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

                    # run cross-validation
                    scores = []
                    for train_idx, val_idx in skf.split(X_train, y_train, groups=None):
                        train_X, val_X = X_train.iloc[train_idx], X_train.iloc[val_idx]
                        train_y, val_y = y_train.iloc[train_idx], y_train.iloc[val_idx]

                        model = lgb.LGBMClassifier(**params, random_state=42)
                        model.fit(train_X, train_y, eval_set=[(val_X, val_y)], early_stopping_rounds=10, verbose=0)
                        
                        y_pred = model.predict(val_X)
                        score = f1_score(val_y, y_pred)
                        scores.append(score)

                    #scores.append(score)
                    
                    # return mean score
                    return np.mean(scores)

                # define hyperparameter ranges
                pbounds = {'n_estimators': (400, 800),
                        'learning_rate': (0.01, 0.1),
                        'max_depth': (3, 9),
                        'num_leaves': (5, 50),
                        'min_child_samples': (10, 100)}

                # run Bayesian optimization
                optimizer = BayesianOptimization(f=lgbm_cv, pbounds=pbounds, random_state=42)
                optimizer.maximize(init_points=10, n_iter=30)

                # print best hyperparameters
                print(optimizer.max)
                f.write('%.3f;' % mean(optimizer.max['target']))
                space = {}
                #scores = cross_val_score(model, atual_XIS, atual_YPSLON[item], scoring='f1', cv=cv, n_jobs=-1)
                for item in optimizer.max['params']:
                    if item != 'learning_rate':
                        space[item] = int(optimizer.max['params'][item])
                    else:
                        space[item] = optimizer.max['params'][item]
                model = lgb.LGBMClassifier(**space, random_state=42)
                model.fit(X_train, y_train)
                #model.fit(X_train, y_train, verbose=0)
                #scores = cross_val_score(model, xis, antibiotic_dfs[antibiotic][antibiotic], scoring='f1', cv=cv, n_jobs=-1)

                y_pred = model.predict(X_test)
                score = f1_score(y_test, y_pred)  
                f.write('%.3f;' % mean(score))
                f.write(str(space))
                # report performance
                #f.write('|'+item +': %.3f (%.3f)' % (mean(scores), std(scores)))
                #f.write('%.3f;' % mean(scores))

            f.write('\n\n')

