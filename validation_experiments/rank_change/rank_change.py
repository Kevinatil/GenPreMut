import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

rbd_name = 'BA.2.1' # 'BA.5.1'

data_root = 'data/rank_change/{}'.format(rbd_name)

# only choose top 100 mutations
sites_dict = {
    "BA.2.1": ['L452R', 'L452Q', 'R346T', 'F486S', 'F486P', 'F486V', 'F490S', 'K444T'],
    "BA.5.1": ['R346T', 'R346S', 'K444T', 'G446R', 'F490S']
}

# the mutation number comparing with wild type RBD
mut_num_dict = {
    "BA.2.1": 16,
    "BA.5.1": 17
}

sites = sites_dict[rbd_name]

df1 = pd.read_csv(os.path.join(data_root, 'mut_rank_all.csv'))
df2 = pd.read_csv(os.path.join(data_root, 'mut_rank_afterexpr.csv'))
df3 = pd.read_csv(os.path.join(data_root, 'mut_rank_final.csv'))

df1 = df1.sort_values(['rank'], ascending = True, axis = 0)['mutation'].values
df2 = df2.sort_values(['rank'], ascending = True, axis = 0)['mutation'].values
df3 = df3.sort_values(['rank'], ascending = True, axis = 0)['mutation'].values

def get_rank(mutations, mut):
    return np.where(mutations == mut)[0][0] + 1

ranks = []

for site in sites:
    ranks.append([get_rank(df1, site), get_rank(df2, site), get_rank(df3, site)])

ranks = np.array(ranks) - mut_num_dict[rbd_name]

print(ranks)


# only visualize top 75 mutation
mask1 = (ranks[:,2] < 75).flatten()
mask2 = (((ranks[:,0] >= ranks[:,2]) | (ranks[:,0] < 5)) & ((ranks[:,1] >= ranks[:,2]) | (ranks[:,2] < 10))).flatten()
mask = mask1 & mask2

ranks = ranks[mask]

plt.figure()
fontsize=15
num = len(ranks)
for i in range(num):
    plt.plot([0,1,2], ranks[i], label = sites[i])
    plt.scatter([0,1,2], ranks[i])
plt.legend(fontsize=fontsize-2)
#plt.legend(fontsize=fontsize-2, ncol=2)
plt.xticks([0,1,2],['stage 1', 'stage 2', 'stage 3'], fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig('ranks_all_{}.svg'.format(rbd_name))
