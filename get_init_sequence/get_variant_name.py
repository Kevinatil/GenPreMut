import os
import re


data_root = "../data"

# lineage_notes.txt is from https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt
f = open(os.path.join(data_root, "raw_seqs/lineage_notes.txt"))

pango2nick = dict()
nick2pango = dict()

for line in f:
    nick, info = line.strip().split('\t')
    if nick == 'Lineage' or '*' in nick:
        continue
    find = re.findall(r'Alias of ([A-Z1-9\.]+)[,\s]', line)
    if len(find):
        pango = find[0]
        pango2nick[pango] = nick
        nick2pango[nick] = pango
    else:
        pango2nick[nick] = nick
        nick2pango[nick] = nick

begin = 'BA.5.1'
end = 'BQ.1'

begin_p = nick2pango[begin]
end_p = nick2pango[end]

print(begin_p, end_p)

pangos = pango2nick.keys()

f = open('candidate_nicks_{}.txt'.format(begin), 'w')
for pango in pangos:
    if pango >= begin_p and pango < end_p:
        f.write('{}\n'.format(pango2nick[pango]))
f.close()