import os
import numpy as np
import pandas as pd

import dmslogo
import dmslogo.colorschemes
import logomaker as lm

os.chdir(os.path.dirname(__file__))


## logo plot code is borrowed from MLAEP(https://github.com/WHan-alter/MLAEP/blob/master/analysis/site_logo_plot.ipynb)


def generate_psedo_count_matrix(pcm_1, pcm_2, len_1, len_2, pseudo_count=0.01):
    # one thing we need to make sure is that the zero counts do not influence our results 
    # As we only focus on the true "mutational" change found by our model
    # so set different pseudo counts according to the size of the two datasets
    pcm_1 = pcm_1 + pseudo_count
    pcm_2 = pcm_2 + pseudo_count*len_2/len_1
    return pcm_1, pcm_2

def get_frequency_matrix(pcm_1, pcm_2):
    # for consistency and avoiding misleading, compute the frequency matrix only from count matrix 
    pfm_1 = pcm_1/np.sum(pcm_1, axis=1)[0]
    pfm_2 = pcm_2/np.sum(pcm_2, axis=1)[0]
    return pfm_1, pfm_2


def kl_logo_matrix(pfm_1, pfm_2, mode="prob_weight_KL"):
    log_pfm_2 = np.log(pfm_2)
    log_pfm_1 = np.log(pfm_1)
    if mode == "prob_weight_KL":
        KL_site = np.sum(pfm_2*np.log(pfm_2) - pfm_2*np.log(pfm_1), axis=1)
        site_norm = np.sum(pfm_2 * np.abs(np.log(pfm_2/pfm_1)), axis=1)[:, np.newaxis]
        kl_logo = KL_site[:, np.newaxis] * pfm_2 * np.log(pfm_2/pfm_1) / site_norm
        #kl_logo = KL_site * pfm_2 * np.log(pfm_2/pfm_1) / site_norm
        return kl_logo
    elif mode=="KL_site":
        KL_site = np.sum(pfm_2*np.log(pfm_2) - pfm_2*np.log(pfm_1), axis=1)
        return KL_site
    


def get_top_df(gisaid_attack,top_n,attack2_init,add):
    gisaid_attack=gisaid_attack.sort_values(by="mean_score",ascending=False)[:top_n]
    counts_mat_init=lm.alignment_to_matrix(attack2_init["seq"])
    counts_mat_attack=lm.alignment_to_matrix(gisaid_attack["seq"])
    counts_mat_init,counts_mat_attack=generate_psedo_count_matrix(counts_mat_init,counts_mat_attack,attack2_init.shape[0],gisaid_attack.shape[0],add)
    # print(counts_mat_init)
    # print(counts_mat_attack)
    pfm_init,pfm_attack2=get_frequency_matrix(counts_mat_init,counts_mat_attack)
    KL=kl_logo_matrix(pfm_init,pfm_attack2,mode="prob_weight_KL")
    KL["pos"]=list(range(331,532))
    df=pd.melt(KL,id_vars="pos")
    df["wt_label"]=df.pos.map(wt_dict)
    df["shade_alpha"]=0.1
    df=df.rename(columns={"value":"prob_weight_KL"})
    return df

def get_KL_site(add,top,gisaid_attack,attack2_init):
    gisaid_attack=gisaid_attack.sort_values(by="mean_score",ascending=False)[:top]
    counts_mat_init=lm.alignment_to_matrix(attack2_init["seq"])
    counts_mat_attack=lm.alignment_to_matrix(gisaid_attack["seq"])
    counts_mat_init,counts_mat_attack=generate_psedo_count_matrix(counts_mat_init,counts_mat_attack,attack2_init.shape[0],gisaid_attack.shape[0],add)
    pfm_init,pfm_attack2=get_frequency_matrix(counts_mat_init,counts_mat_attack)
    KL_site=kl_logo_matrix(pfm_init,pfm_attack2,"KL_site")
    KL_site.index=KL_site.index+331
    return(KL_site)


name = 'BA.2.1'
#name = 'BA.5.1'

data_root = 'data/'

# drawing without 346 for subplot
with346 = True

if with346:
    targets_all_dict = {
        'BA.2.1': ['452R','452Q','346T','486S','486P','460K','486V','446S','490S','339H','444T','445P','368I','478R'],
        'BA.5.1': ['346T','346S','444T','446R','460K','490S']
    }

    # rank top 100
    targets_todraw_dict = {
        'BA.2.1': ['452R','452Q','346T','446S','490S'],
        'BA.5.1': ['346T','460K']
    }
else:
    targets_all_dict = {
        'BA.2.1': ['452R','452Q','486S','486P','460K','486V','446S','490S','339H','444T','445P','368I','478R'],
        'BA.5.1': ['444T','446R','460K','490S']
    }

    # rank top 100
    targets_todraw_dict = {
        'BA.2.1': ['452R','452Q','446S','490S'],
        'BA.5.1': ['460K']
    }


synthetic_init=pd.read_csv(os.path.join(data_root, 'logo/MLAEP/origin_{}.csv'.format(name)), index_col=0)
synthetic=pd.read_csv(os.path.join(data_root, 'logo/MLAEP/synthesize_{}.csv'.format(name)), index_col=0)

with open(os.path.join(data_root, 'logo/MLAEP/Covid19_RBD_seq.txt'),'r') as f:
    wt_seq=f.read()
wt=[str(x)+str(y) for x,y in zip(list(wt_seq),list(range(331,532)))]
wt_dict=pd.DataFrame({"pos":list(range(331,532)),"wt":wt}).set_index("pos").to_dict()
wt_dict=wt_dict["wt"]

## set the pseudo counts and the sequence number selected
add=0.1
top=synthetic.shape[0]

KL_site=get_KL_site(add,top,synthetic,synthetic_init)
KL_site=KL_site.to_frame(name="KL_site")

## get the dataframe contain all the information
df=get_top_df(synthetic,top,synthetic_init,add)
df=df.merge(KL_site,on="pos")

df['color'] = '#808080'


targets_all = targets_all_dict[name]
targets_todraw = targets_todraw_dict[name]

sites_todraw = []
for site in targets_all:
    sites_todraw.append(int(site[:3]))


for site in targets_todraw:
    pos,mut = int(site[:3]), site[3]
    idx = df[(df['pos'] == pos)&(df['variable'] == mut)].index
    assert len(idx)
    df.loc[idx,'color'] = '#610345'


mask = []
for i in range(len(df)):
    mask.append(df.loc[i]['wt_label'][0] != df.loc[i]['variable'])
mask = np.array(mask)
df = df[mask]

df = df.sort_values(by = ['pos','prob_weight_KL'], ascending=False)
df.index = np.arange(len(df))

print(df.head())

df['prob_weight_KL'] = df['prob_weight_KL']-df['prob_weight_KL'].min()/(df['prob_weight_KL'].max()-df['prob_weight_KL'].min())
df['KL_site'] = df['KL_site']-df['KL_site'].min()/(df['KL_site'].max()-df['KL_site'].min())


idxs = []
for i in range(331, 532):
    df_ = df[df['pos']==i][:10].index
    idxs += df_.tolist()

df = df.loc[idxs]

fig, ax = dmslogo.draw_logo(data=df[df["pos"].isin(sites_todraw)],
                        x_col='pos',
                        letter_col='variable',
                        letter_height_col="prob_weight_KL",
                        color_col="color",
                        xtick_col="wt_label",
                        xlabel="Site",
                        ylabel="Bits",
                        axisfontscale= 3,
                        letterheightscale=1,
                )


if with346:
    fig.savefig('{}_KL.svg'.format(name),bbox_inches='tight')
else:
    fig.savefig('{}_KL_no346.svg'.format(name),bbox_inches='tight')