from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import  *
from sklearn.model_selection import GridSearchCV
import pickle
from aaindex import aaindex1
from tqdm import tqdm
from FFT import get_aai_encoding
import numpy as np
import pandas as pd

def regscore(fftdata,target,cv,return_index_model=False):
    X = fftdata
    y= target
    reg = PLSRegression()
    param_grid = {'n_components': range(2, 11)}


    gsearch = GridSearchCV(reg, param_grid,cv=cv,n_jobs=-1,scoring='neg_mean_squared_error')
    try:
        gsearch.fit(X, y)
        best_model=gsearch.best_estimator_
        y_pred = gsearch.predict(X)
        mse = mean_squared_error( y,y_pred)
        r2 = r2_score( y,y_pred)
    except ValueError as e:  
        mse=100
        r2=-100
        gsearch.best_score_=-100
        gsearch.best_params_={'n_components':1}

    
    if return_index_model:
        pickle.dump(best_model,open(f'./{return_index_model}.dat','wb'))

    return mse,r2,-gsearch.best_score_,gsearch.best_params_['n_components']


def index_search(mutdict,target,indexlist,cv):
    score ={}
    all_indices = aaindex1.record_codes()

    progress_bar = tqdm(all_indices, desc='Processing')
    for indexname in progress_bar:
        indexlist_=indexlist.copy()
        indexlist_.append(indexname)
        fftdata= get_aai_encoding(mutdict,indices=indexlist_) 
        indexname_ = '_'.join(indexlist_)
        score[indexname_] = regscore(fftdata,target,cv=cv)
    score_sort = pd.DataFrame(score,index=['MSE','R2','cvMSE','n_components']).transpose().sort_values(by='cvMSE', ascending=True)
    score_sort['n_components'] = score_sort['n_components'].astype(int)
    screen_index =score_sort.index.values[0]
    return screen_index.split('_'),score_sort

def make_prediction(train_dict,target,predict_dict,indexlist,n_idx): 
    indexlist_ =indexlist.copy()
    indexname= '_'.join(indexlist_)
    train = get_aai_encoding(train_dict,indices=indexlist_)
    predict = get_aai_encoding(predict_dict,indices=indexlist_)
    print('The number of mutants for training:',len(train))
    print('The number of mutants for prediction:',len(predict))
    X=train
    y=target
    reg = PLSRegression(n_components=n_idx)
    all_pred = reg.fit(X,y).predict(predict)
    result=pd.DataFrame(index =[key for key in predict_dict.keys()]) 
    result[indexname]=all_pred
    return result