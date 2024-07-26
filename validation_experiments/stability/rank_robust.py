import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

data_root = 'data/'

# calculate the average rankings
if 0:
    # three repeat experiments
    idx = 0 # 1, 2

    df_old=pd.read_csv(os.path.join(data_root, 'repeat_experiments/mut_freq_model1_{}.csv'.format(idx)))
    df_new=pd.read_csv(os.path.join(data_root, 'repeat_experiments/mut_freq_model2_{}.csv'.format(idx)))

    muts_old=set(df_old['mutation'].values)
    muts_new=set(df_new['mutation'].values)

    df_old = df_old['mutation'].values
    df_new = df_new['mutation'].values

    muts = muts_old.intersection(muts_new)

    rank_old = []
    rank_new = []

    for i in range(len(df_old)):
        if df_old[i] in muts:
            rank_old.append(df_old[i])

    for i in range(len(df_new)):
        if df_new[i] in muts:
            rank_new.append(df_new[i])

    rank_old = np.array(rank_old)
    rank_new = np.array(rank_new)

    rank_all = []

    for mut in muts:
        rank_all.append([mut, (np.where(rank_old == mut)[0][0] + np.where(rank_new == mut)[0][0]) / 2])

    rank_all.sort(key = lambda k: k[1])

    pd.DataFrame(rank_all, columns = ['mutation', 'rank']).to_csv(os.path.join(data_root, 'repeat_experiments/mut_rank{}.csv'.format(idx)), index=False)


# BA.2.1 as parent node
ori_muts = ['G339D','S371F','S373P','S375F','T376A','D405N','R408S','K417N','N440K','S477N','T478K','E484A','Q493R','Q498R','N501Y','Y505H']



df0 = pd.read_csv(os.path.join(data_root, 'repeat_experiments/mut_rank0.csv'))['mutation'].values
df1 = pd.read_csv(os.path.join(data_root, 'repeat_experiments/mut_rank1.csv'))['mutation'].values
df2 = pd.read_csv(os.path.join(data_root, 'repeat_experiments/mut_rank2.csv'))['mutation'].values

fontsize = 20


# 0 1
rank01 = []
for mut in df0:
    if mut in ori_muts:
        continue
    tp_ = np.where(df1==mut)[0]
    if len(tp_):
        rank01.append([np.where(df0==mut)[0][0], tp_[0]])

rank01 = np.array(rank01)
print(np.corrcoef(rank01[:,0], rank01[:,1]))

if 1:
    rank01 = np.array(rank01)
    plt.figure(figsize=(10, 8))
    plt.scatter(rank01[:,1], rank01[:,0], s=5)
    plt.plot([10,np.max(rank01)],[10,np.max(rank01)],linestyle='--',linewidth=2,color='orange')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig('rank01.svg')

    plt.figure(figsize=(10, 8))
    plt.scatter(rank01[:,1], rank01[:,0], s=5)
    plt.plot([10,np.max(rank01)],[10,np.max(rank01)],linestyle='--',linewidth=2,color='orange')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(left = -10, right=300)
    plt.ylim(bottom = -10, top=300)
    #plt.scatter(100,100,s=20,c='red')
    #plt.text(100,80,'rank 100',fontsize=fontsize)
    #plt.plot([10,10],[10,180],linestyle='--',linewidth=2,color='red')
    #plt.plot([100,100],[100,180],linestyle='--',linewidth=2,color='red')
    plt.savefig('rank01_300.svg')


# 0 2
rank02 = []
for mut in df0:
    if mut in ori_muts:
        continue
    tp_ = np.where(df2==mut)[0]
    if len(tp_):
        rank02.append([np.where(df0==mut)[0][0], tp_[0]])

rank02 = np.array(rank02)
print(np.corrcoef(rank02[:,0], rank02[:,1]))


if 1:
    rank02 = np.array(rank02)
    plt.figure(figsize=(10, 8))
    plt.scatter(rank02[:,1], rank02[:,0], s=5)
    plt.plot([10,np.max(rank02)],[10,np.max(rank02)],linestyle='--',linewidth=2,color='orange')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig('rank02.svg')

    plt.figure(figsize=(10, 8))
    plt.scatter(rank02[:,1], rank02[:,0], s=5)
    plt.plot([10,np.max(rank02)],[10,np.max(rank02)],linestyle='--',linewidth=2,color='orange')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(left = -10, right=300)
    plt.ylim(bottom = -10, top=300)
    #plt.scatter(100,100,s=20,c='red')
    #plt.text(100,80,'rank 100',fontsize=fontsize)
    #plt.plot([10,10],[10,180],linestyle='--',linewidth=2,color='red')
    #plt.plot([100,100],[100,180],linestyle='--',linewidth=2,color='red')
    plt.savefig('rank02_300.svg')

