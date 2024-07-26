# data process for the training and prediction of bind, expr
import os
import numpy as np

import pandas as pd

import random

data_root = 'data'

# bind data
if 1:
    df_bind=pd.read_csv(os.path.join(data_root, 'property_prediction/bind/bind_variant_data.csv'))
    # only use samples of BA.1 and BA.2
    mask_bind=(df_bind['target']=='Omicron_BA1')|(df_bind['target']=='Omicron_BA2')
    label_bind=df_bind[mask_bind]['bind_ratio'].values
    # load features extracted by backbone
    feature_bind=np.load(os.path.join(data_root, 'property_prediction/embedding_all_data.npy'))
    feature_bind=np.array(feature_bind)

    assert len(feature_bind)==len(label_bind)

    mask=np.isnan(label_bind)
    feature_bind=feature_bind[~mask]
    label_bind=label_bind[~mask]

    print(feature_bind.shape)
    print(label_bind.shape)
    np.save(os.path.join(data_root, 'property_prediction/bind/all_bind_data.npy'),feature_bind)
    np.save(os.path.join(data_root, 'property_prediction/bind/all_bind_label.npy'),label_bind)

    del feature_bind

# expr data
if 1:
    df_expr=pd.read_csv(os.path.join(data_root, 'property_prediction/expr/expr_variant_data.csv'))
    # only use samples of BA.1 and BA.2
    mask_expr=(df_expr['target']=='Omicron_BA1')|(df_expr['target']=='Omicron_BA2')
    label_expr=df_expr[mask_expr]['expr_ratio'].values
    # load features extracted by customized backbone
    feature_expr=np.load(os.path.join(data_root, 'property_prediction/embedding_all_data.npy'))
    feature_expr=np.array(feature_expr)

    assert len(feature_expr)==len(label_expr)

    mask=np.isnan(label_expr)
    feature_expr=feature_expr[~mask]
    label_expr=label_expr[~mask]

    print(feature_expr.shape)
    print(label_expr.shape)
    np.save(os.path.join(data_root, 'property_prediction/expr/all_expr_data.npy'),feature_expr)
    np.save(os.path.join(data_root, 'property_prediction/expr/all_expr_label.npy'),label_expr)

###################################################################


###################################################################
# train test data split

# We use the same data splitting as in E2VD, please download the splitting files as https://github.com/ZhiweiNiepku/E2VD.
# For reproduction, please download the extracted features and splitting masks to data_root, or download the trained checkpoints from E2VD directly.
if 0:
    #bind
    label=np.load(os.path.join(data_root, 'property_prediction/bind/all_bind_label.npy'))
    mask_=np.where(np.isnan(label))[0]
    pidx=random.sample(np.where(label>=1)[0].tolist(),79)
    nidx=random.sample(np.where(label<1)[0].tolist(),79)
    test_idx=np.array(pidx+nidx)
    train_idx=[]
    for i in range(len(label)):
        if i not in test_idx and i not in mask_:
            train_idx.append(i)
    train_idx=np.array(train_idx)

    print(train_idx.shape)
    print(test_idx.shape)
    np.save(os.path.join(data_root, 'property_prediction/bind/bind_train_idx.npy'),train_idx)
    np.save(os.path.join(data_root, 'property_prediction/bind/bind_test_idx.npy'),test_idx)

    #expr
    label=np.load(os.path.join(data_root, 'property_prediction/expr/all_expr_label.npy'))
    mask_=np.where(np.isnan(label))[0]
    pidx=random.sample(np.where(label>=1)[0].tolist(),90)
    nidx=random.sample(np.where(label<1)[0].tolist(),90)
    test_idx=np.array(pidx+nidx)
    train_idx=[]
    for i in range(len(label)):
        if i not in test_idx and i not in mask_:
            train_idx.append(i)
    train_idx=np.array(train_idx)

    print(train_idx.shape)
    print(test_idx.shape)
    np.save(os.path.join(data_root, 'property_prediction/expr/expr_train_idx.npy'),train_idx)
    np.save(os.path.join(data_root, 'property_prediction/expr/expr_test_idx.npy'),test_idx)


if 1:
    feature_all=np.load(os.path.join(data_root, 'property_prediction/bind/all_bind_data.npy'))
    label_all=np.load(os.path.join(data_root, 'property_prediction/bind/all_bind_label.npy'))
    print(feature_all.shape)
    train_idx=np.load(os.path.join(data_root, 'property_prediction/bind/bind_train_idx.npy'))
    test_idx=np.load(os.path.join(data_root, 'property_prediction/bind/bind_test_idx.npy'))
    train_feature=feature_all[train_idx]
    train_label=label_all[train_idx]

    test_feature=feature_all[test_idx]
    test_label=label_all[test_idx]

    print(train_feature.shape)
    print(train_label.shape)
    print(test_feature.shape)
    print(test_label.shape)

    np.save(os.path.join(data_root, 'property_prediction/bind/bind_train_data.npy'), train_feature)
    np.save(os.path.join(data_root, 'property_prediction/bind/bind_train_label.npy'), train_label)
    np.save(os.path.join(data_root, 'property_prediction/bind/bind_test_data.npy'), test_feature)
    np.save(os.path.join(data_root, 'property_prediction/bind/bind_test_label.npy'), test_label)



    feature_all=np.load(os.path.join(data_root, 'property_prediction/expr/all_expr_data.npy'))
    label_all=np.load(os.path.join(data_root, 'property_prediction/expr/all_expr_label.npy'))
    print(feature_all.shape)
    train_idx=np.load(os.path.join(data_root, 'property_prediction/expr/expr_train_idx.npy'))
    test_idx=np.load(os.path.join(data_root, 'property_prediction/expr/expr_test_idx.npy'))
    train_feature=feature_all[train_idx]
    train_label=label_all[train_idx]

    test_feature=feature_all[test_idx]
    test_label=label_all[test_idx]

    print(train_feature.shape)
    print(train_label.shape)
    print(test_feature.shape)
    print(test_label.shape)

    np.save(os.path.join(data_root, 'property_prediction/expr/expr_train_data.npy'), train_feature)
    np.save(os.path.join(data_root, 'property_prediction/expr/expr_train_label.npy'), train_label)
    np.save(os.path.join(data_root, 'property_prediction/expr/expr_test_data.npy'), test_feature)
    np.save(os.path.join(data_root, 'property_prediction/expr/expr_test_label.npy'), test_label)
