import os
import numpy as np
import pandas as pd

import dmslogo
import dmslogo.colorschemes

name = 'BA.2.1'
# name = 'BA.5.1'

data_root = 'data/'

targets_all_dict = {
    'BA.2.1': ['452R','452Q','346T','486S','486P','460K','486V','446S','490S','339H','444T','445P','368I','478R'],
    'BA.5.1': ['346T','346S','444T','446R','460K','490S']
}

# rank top 100
targets_todraw_dict = {
    'BA.2.1': ['452R','452Q','346T','486S','486P','460K','486V','446S','490S','339H','444T'],
    'BA.5.1': ['346T','346S','444T','446R','460K','490S']
}

df=pd.read_csv(os.path.join(data_root, 'logo/GenPreMut/logo_draw_{}.csv'.format(name)))

df['pos']=df['mutation'].apply(lambda x: x[1:4]).apply(int)
df['variable']=df['mutation'].apply(lambda x: x[4])
df['wt_label']=df['mutation'].apply(lambda x: x[:4])
df['color']='#808080'

freq = df.groupby('pos').sum()['frequency']
df['freq_sum'] = df['pos'].apply(lambda x: freq[x])


targets_all = targets_all_dict[name]
targets_todraw = targets_todraw_dict[name]


sites_todraw = []
for site in targets_all:
    sites_todraw.append(int(site[:3]))

for site in targets_todraw:
    pos,mut = site[:3], site[3]
    idx = df[(df['pos'] == int(pos))&(df['variable'] == mut)].index
    assert len(idx)
    df.loc[idx,'color'] = '#610345'



fig, ax = dmslogo.draw_logo(data=df[df["pos"].isin(sites_todraw)],
                        x_col='pos',
                        letter_col='variable',
                        letter_height_col='frequency',
                        color_col="color",
                        xtick_col="wt_label",
                        xlabel="Site",
                        ylabel="frequency",
                        axisfontscale= 2,
                        letterheightscale=1,
                )

fig.savefig('{}_freq.svg'.format(name), bbox_inches='tight')