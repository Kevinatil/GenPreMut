import os
import numpy as np
import pandas as pd

from tqdm import tqdm

rbd_name = 'BA.2.1'
data_root = 'data/'

wt='NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST'

def get_mut(seq):
    ret=[]
    for i in range(len(wt)):
        if wt[i]!=seq[i]:
            ret.append('{}{}{}'.format(wt[i],i+331,seq[i]))
    return ", ".join(ret)

filen=100

df=[]

for i in tqdm(range(filen)):
    f = open(os.path.join(data_root, 'raw_seqs/mutation_{}_{}.txt'.format(rbd_name, i)))
    bind_cls = np.load(os.path.join(data_root, 'predicts/bind/{}/predict_cls_all_{}.npy'.format(rbd_name,i)))
    expr_cls = np.load(os.path.join(data_root, 'predicts/expr/{}/predict_cls_all_{}.npy'.format(rbd_name,i)))
    
    lines_=f.readlines()
    f.close()
    
    for j in range(len(lines_)):
        df.append([lines_[j].strip(), get_mut(lines_[j].strip()), bind_cls[j], expr_cls[j]])

df=pd.DataFrame(df,columns=['sequence', 'mutation', 'bind_cls', 'expr_cls'])

print(df.shape)
print(df.head())
df.to_csv(os.path.join(data_root, 'predicts/df/pred_{}.csv'.format(rbd_name)),index=False)