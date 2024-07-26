import os
import numpy as np
import pandas as pd

import pickle
from tqdm import tqdm

# model1 use_res_clean.csv, supplementary2.csv download: https://github.com/jianfcpku/convergent_RBD_evolution

model_path = 'ckpt'

rbd_name = 'BA.2.1' # 'BA.5.1'

if rbd_name == 'BA.2.1':
    ab_use = ['BA.1 convalescents', 'BA.2 convalescents']
else:
    ab_use = ['BA.1 convalescents', 'BA.2 convalescents', 'BA.5 convalescents']

# choose antibodies to use
if 1:
    df=pd.read_csv('supplementary2.csv',skiprows=[0])
    mask_ = np.zeros(len(df)).astype(bool)
    for source in ab_use:
        mask_ |= (df['source']==source).values
    antibody_names=np.unique(df[mask_]['Antibody  Name'])

    print(len(antibody_names))

    df=pd.read_csv('use_res_clean.csv')
    antibodies=df['antibody'].values
    mask_=[]
    for i in range(len(antibodies)):
        mask_.append((antibodies[i] in antibody_names))

    df=df[mask_]

    print(df.shape)

    df.to_csv(os.path.join(model_path, 'antibody_barrier_model/model1/escape_antibody_{}.csv'.format(rbd_name)), index=False)




df=pd.read_csv(os.path.join(model_path, 'antibody_barrier_model/model1/escape_antibody_{}.csv'.format(rbd_name)))
groups=df['group'].unique()
antibody_groups=df[['antibody','group']]
groupby=antibody_groups.groupby('group')

group_dict={}

for i in range(len(groups)):
    group_dict[groups[i]]=groupby.get_group(groups[i])['antibody'].unique()



# get group_mean.pkl
if 1:
    res_dict={}

    groupby=df.groupby('group')
    for i in range(len(groups)):
        
        group_=groupby.get_group(groups[i])
        sites=group_['site'].unique()
        groupby_site=group_.groupby(['site'])

        res_site_dict={}
        for s in sites:
            
            group_site_=groupby_site.get_group(s)
            muts=group_site_['mutation'].unique()
            group_muts=group_site_.groupby(['mutation'])

            res_mut_dict={}
            for m in muts:
                groupby_muts_=group_muts.get_group(m)
                res_mut_dict[m]=groupby_muts_['mut_escape'].sum()/len(group_dict[groups[i]])
            res_site_dict[s]=res_mut_dict
        res_dict[groups[i]]=res_site_dict

    with open(os.path.join(model_path, 'antibody_barrier_model/model1/group_mean_{}.pkl'.format(rbd_name)), 'wb') as f:
        pickle.dump(res_dict,f)
