import os
import pandas as pd


name = 'BA.2.1' # 'BA.5.1'
data_root = 'data/'

df = pd.read_csv(os.path.join(data_root, 'logo/MLAEP/success_{}_seq.txt'.format(name)), 
                 index_col=0, sep='\t', 
                 names=['seq_ori','col1','col2','seq','col3','mean_score'])

df[['seq', 'mean_score']].to_csv(os.path.join(data_root, 'logo/MLAEP/synthesize_{}.csv'.format(name)), index=False)

