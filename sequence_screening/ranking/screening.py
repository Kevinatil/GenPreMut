import os
import re
import pandas as pd


rbd_name = 'BA.2.1'
data_root = 'data/'

ori='NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST'


def get_muts(seq):
    len_=len(seq)
    ret=[]
    for i in range(len_):
        if ori[i]!=seq[i]:
            ret.append('{}{}{}'.format(ori[i],331+i,seq[i]))
    return ret

for kind in ['model1','model2']:
    thres=0.5
    df=pd.read_csv(os.path.join(data_root, 'predicts/df/pred_{}.csv'.format(rbd_name)))

    len_=len(df)

    seq1=set(df[df['expr_cls']>thres]['sequence'].values)
    seq2=set(df.sort_values(['antibody_barrier_{}'.format(kind)],ascending=False)['sequence'].values[:len_//2])
    seqs=seq1.intersection(seq2)
    seqs=list(seqs)

    f=open(os.path.join(data_root, 'predicts/df/seqs_final_{}_{}.txt'.format(rbd_name, kind)), 'w')
    for seq in seqs:
        f.write('{}\n'.format(seq))
    f.close()
    print(len(seqs))

    all_muts=[]
    for seq in seqs:
        all_muts+=get_muts(seq)

    mut_types=set(all_muts)
    nums={}
    for i in mut_types:
        nums[i]=0
    for i in all_muts:
        nums[i]+=1

    freqs=[]
    for i in nums.keys():
        freqs.append([i,nums[i]])
    df=pd.DataFrame(freqs,columns=['mutation','frequency'])
    df.sort_values(['frequency'],ascending=False).to_csv(os.path.join(data_root, 'predicts/df/mut_freq_{}_{}.csv'.format(rbd_name, kind)), index=False)