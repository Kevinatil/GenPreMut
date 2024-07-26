import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df1 = pd.DataFrame([[1000000, 373201, 213784], [1000000, 373313, 214045], [1000000, 373418, 213597]], columns = ['all', 'after expression', 'after barrier I'], index=['experiment 1', 'experiment 2', 'experiment 3']).T

df2 = pd.DataFrame([[1000000, 373201, 179474], [1000000, 373313, 179519], [1000000, 373418, 179477]], columns = ['all', 'after expression', 'after barrier II'], index=['experiment 1', 'experiment 2', 'experiment 3']).T

ticks = ['0', '0.2M', '0.4M', '0.6M', '0.8M', '1M']

fontsize = 15

if 0:
    plt.figure(figsize=(8, 8))
    df1.plot.bar()
    plt.xticks(rotation=0, fontsize = fontsize)
    plt.yticks(np.linspace(0, 1, len(ticks)) * 1e6, ticks, fontsize = fontsize)
    plt.legend(fontsize=fontsize)
    plt.savefig('bar1.svg')


    plt.figure(figsize=(8, 8))
    df2.plot.bar()
    plt.xticks(rotation=0, fontsize = fontsize)
    plt.yticks(np.linspace(0, 1, len(ticks)) * 1e6, ticks, fontsize = fontsize)
    plt.legend(fontsize=fontsize)
    plt.savefig('bar2.svg')
    
if 1:
    df1 = pd.DataFrame([[1000000, 373201, 213784], [1000000, 373313, 214045], [1000000, 373418, 213597]], columns = ['all', 'after expression', 'after barrier I'], index=['exp 1', 'exp 2', 'exp 3'])
    df1['all'] = df1['all'] - df1['after expression']
    df1['after expression'] = df1['after expression'] - df1['after barrier I']
    df1 = df1[['after barrier I', 'after expression', 'all']]


    
    plt.figure()
    df1.plot.barh(#figsize=(5, 4), 
                  stacked=True, color={"all": "#95D0BF", "after expression": "#00A5B9", "after barrier I": "#055C8D"})
    plt.yticks(rotation=90, fontsize = fontsize)
    plt.xticks(np.linspace(0, 1, len(ticks)) * 1e6, ticks, fontsize = fontsize)
    plt.ylim(top=2.8)
    plt.legend(fontsize=fontsize-3.5,ncol=3)
    plt.savefig('bar1s.svg')

    
    df2 = pd.DataFrame([[1000000, 373201, 179474], [1000000, 373313, 179519], [1000000, 373418, 179477]], columns = ['all', 'after expression', 'after barrier II'], index=['exp 1', 'exp 2', 'exp 3'])
    df2['all'] = df2['all'] - df2['after expression']
    df2['after expression'] = df2['after expression'] - df2['after barrier II']
    df2 = df2[['after barrier II', 'after expression', 'all']]

    plt.figure()
    df2.plot.barh(#figsize=(5, 4), 
                  stacked=True, color={"all": "#95D0BF", "after expression": "#00A5B9", "after barrier II": "#055C8D"})
    plt.yticks(rotation=90, fontsize = fontsize)
    plt.xticks(np.linspace(0, 1, len(ticks)) * 1e6, ticks, fontsize = fontsize)
    plt.ylim(top=2.8)
    plt.legend(fontsize=fontsize-3.5,ncol=3)
    plt.savefig('bar2s.svg')