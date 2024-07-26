import os
import numpy as np
import pandas as pd

import pickle
from tqdm import tqdm


model_root = 'ckpt/antibody_barrier_model'

def get_escape_score(group_mean, seq, group):
    len_ = len(seq)
    ret = 0
    for i in range(len_):
        try:
            s = group_mean[group][i+331][seq[i]]
        except:
            s = 0
        ret += s
    return ret

def get_group_weight(group_dict, counter, key):
    return len(group_dict[key])/counter

def antibody_barrier(df_path, rbd_name, version):
    root = os.path.join(model_root, 'model{}'.format(version))

    df_ab = pd.read_csv(os.path.join(root, 'escape_antibody_{}.csv'.format(rbd_name)))
    groups = df_ab['group'].unique().tolist()
    antibody_groups = df_ab[['antibody','group']]
    groupby = antibody_groups.groupby('group')

    group_dict = dict()
    for i in range(len(groups)):
        group_dict[groups[i]] = groupby.get_group(groups[i])['antibody'].unique()

    group_mean = pickle.load(open(os.path.join(root, 'group_mean_{}.pkl'.format(rbd_name)),'rb'))


    df=pd.read_csv(df_path)
    seqs=df['sequence'].values
    res=[]
    for seq in tqdm(seqs):
        res_=[]
        for group_ in groups:
            res_.append(get_escape_score(group_mean, seq, group_))
        res.append(res_)

    df = pd.DataFrame()
    df[groups] = res

    counter=0
    for key in group_dict.keys():
        counter+=len(group_dict[key])
    print(counter)

    scores=[]
    for i in tqdm(range(len(df))):
        score_=0
        for group in groups:
            score_ += get_group_weight(group_dict, counter, group) * df[group][i]
        scores.append(score_)

    return scores


if __name__ == "__main__":

    rbd_name = 'BA.2.1'

    df_path = 'data/predicts/df/pred_{}.csv'.format(rbd_name)
    df = pd.read_csv(df_path)

    df['antibody_barrier_score1'] = antibody_barrier(df_path, rbd_name, version=1)
    df['antibody_barrier_score2'] = antibody_barrier(df_path, rbd_name, version=2)

    df.to_csv(df_path, index=False)