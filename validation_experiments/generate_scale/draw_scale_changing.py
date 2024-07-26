import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from tqdm import tqdm

data_root = 'data'

files=['50K','100K','150K','200K','250K','300K','350K','400K','450K',
       '500K','550K','600K','650K','700K','750K','800K','850K','900K',
       '950K','1M','1.050M','1.1M','1.15M','1.2M']
process=24
per_num=50000

df_m1_all=[]
df_m1=pd.read_csv(os.path.join(data_root, 'generate_scale/full_qab_model1.csv'))
for i in range(process):
    tp_df=df_m1[:(i+1)*per_num]
    df_m1_all.append(tp_df)


df_m2_all=[]
df_m2=pd.read_csv(os.path.join(data_root, 'generate_scale/full_qab_model2.csv'))
for i in range(process):
    tp_df=df_m2[:(i+1)*per_num]
    df_m2_all.append(tp_df)



nums=[5000,10000,15000,20000]

means_mean_m1_all=[]
for i in range(len(files)):
    tp_df=df_m1_all[i].sort_values(['mean'], ascending=False)
    
    means_mean=[]
    for j in nums:
        means_mean.append(tp_df['mean'][:j].mean())
    means_mean_m1_all.append(means_mean)

    #print(means_mean)

means_mean_m1_all = np.array(means_mean_m1_all).T




means_mean_m2_all=[]
for i in range(len(files)):
    tp_df=df_m2_all[i].sort_values(['mean'], ascending=False)
    
    means_mean=[]
    for j in nums:
        means_mean.append(tp_df['mean'][:j].mean())
    means_mean_m2_all.append(means_mean)

    #print(means_mean)

means_mean_m2_all = np.array(means_mean_m2_all).T

# merge
means_mean_m1_all_dict=dict()
for i in range(len(nums)):
    means_mean_m1_all_dict[nums[i]]=means_mean_m1_all[i]
    
means_mean_m2_all_dict=dict()
for i in range(len(nums)):
    means_mean_m2_all_dict[nums[i]]=means_mean_m2_all[i]


plt.figure(figsize=(10,8))
num_to_draw=nums
for i in range(len(num_to_draw)):
    means_mean=means_mean_m1_all_dict[num_to_draw[i]]
    deltas=[]
    for j in range(len(means_mean)-1):
        deltas.append(means_mean[j+1]-means_mean[j])
    x_ = np.arange(len(deltas)) + 0.5
    plt.scatter(x_, deltas) #,label='{}'.format(nums[i]))
    plt.plot(x_, deltas,label='top {}'.format(num_to_draw[i]))

fontsize = 20
plt.xticks(np.arange(len(files))[::4],files[::4],fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.plot([0,23.2],[0,0],'r--')

ymin,ymax=plt.ylim()
plt.ylim(ymin,0.029)

xmin,xmax=plt.xlim()
plt.xlim(xmin,23.5)

plt.legend(fontsize=fontsize)
plt.savefig('qab_change_m1.svg')



# new
plt.figure(figsize=(10,8))
num_to_draw=nums
for i in range(len(num_to_draw)):
    means_mean=means_mean_m2_all_dict[num_to_draw[i]]
    deltas=[]
    for j in range(len(means_mean)-1):
        deltas.append(means_mean[j+1]-means_mean[j])
    plt.scatter(np.arange(len(deltas)),deltas) #,label='{}'.format(nums[i]))
    plt.plot(deltas,label='top {}'.format(num_to_draw[i]))

plt.xticks(np.arange(len(files))[::4],files[::4],fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.plot([0,23.2],[0,0],'r--')

xmin,xmax=plt.xlim()
plt.xlim(xmin,23.5)

plt.legend(fontsize=fontsize)
plt.savefig('qab_change_m2.svg')
