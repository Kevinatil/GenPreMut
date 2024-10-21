import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

name = 'BA.2.1'
#name = 'BA.5.1'

# rank top 100
targets_todraw_dict = {
    'BA21': ['452R','452Q','346T','486S','486P','460K','486V','446S','490S','339H','444T'],
    'BA51': ['346T','346S','444T','446R','460K','490S','450D']
}

if name == 'BA21':
    ranks_ms = {
        '452R':1,
        '452Q':3,
        '346T':8,
        '486S':9,
        '486P':13,
        '460K':27,
        '486V':32,
        '446S':35,
        '490S':40,
        '339H':62,
        '444T':73
    }

    ranks_mlaep = {
        '452R':59,
        '452Q':73,
        '346T':22,
        '486S':2397,
        '486P':2768,
        '460K':1366,
        '486V':2376,
        '446S':71,
        '490S':76,
        '339H':218,
        '444T':851
    }
else:
    ranks_ms = {
        '346T':2,
        '346S':3,
        '444T':35,
        '446R':23,
        '460K':31,
        '490S':22,
        '450D':74
    }

    ranks_mlaep = {
        '346T':2,
        '346S':3606,
        '444T':2578,
        '446R':3364,
        '460K':276,
        '490S':3615,
        '450D':341
    }



if name == 'BA.2.1':
    plt.figure()
    fontsize=15
    size = 50
if name == 'BA.5.1':
    plt.figure(figsize=(6, 8))
    fontsize=25
    size = 120
sites = targets_todraw_dict[name]
rank_ms = []
rank_mlaep = []
for site in sites:
    num_ = ranks_ms[site] if ranks_ms[site] <= 100 else 120
    rank_ms.append(num_)
    num_ = ranks_mlaep[site] if ranks_mlaep[site] <= 100 else 120
    rank_mlaep.append(num_)

plt.scatter(np.arange(len(rank_ms)), rank_ms, s=size, label='MutSeek', zorder=2)
plt.scatter(np.arange(len(rank_mlaep)), rank_mlaep, s=size, label='MLAEP', zorder=2)

for i in range(len(sites)):
    plt.plot([i, i], [rank_ms[i], rank_mlaep[i]], color = 'grey', linewidth=5, zorder=1)

plt.xticks(np.arange(len(sites)), sites, fontsize = fontsize, rotation=90)
plt.yticks(np.arange(0, 140, 20), ['0','20','40','60','80','100','>100'], fontsize = fontsize)

plt.legend(fontsize=fontsize)

plt.savefig('rank_change_{}.svg'.format(name))
