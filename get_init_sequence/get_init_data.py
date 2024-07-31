import os
import re
import glob

import time
import numpy as np
from tqdm import tqdm

from Bio import SeqIO

data_root = "data"
model_root = "ckpt"

rbd_name = "BA.2.1" # "BA.5.1" # "XBB.1.5"

parent_node = {
    "BA.2.1": "NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "BA.5.1": "NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVGGNYNYRYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGVNCYFPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "XBB.1.5":"NITNLCPFHEVFNATTFASVYAWNRKRISNCVADYSVIYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKPSGNYNYLYRLFRKSKLKPFERDISTEIYQAGNKPCNGVAGPNCYSPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST"
}

## filter csv to use
csv_path=os.path.join(data_root, 'GISAID_ID', 'total_csv')
all_files=os.listdir(csv_path)

def file_filter(name, nicks, rbd_name):
    if rbd_name == "XBB.1.5":
        return "XBB" in name
    else:
        return name in nicks

f = open('candidate_nicks_{}.txt'.format(rbd_name), 'r')
nicks = f.readlines()
f.close()
nicks = [nick.strip() for nick in nicks]

filtered_file=[]
for file in all_files:
    if file_filter(file, nicks, rbd_name):
        filtered_file.append(file)
        
        
if 0:
    # The downloaded fasta file may contain bad lines, which will cause parsing error. 
    # If so, please delete bad lines before parsing.
    
    # Please download the latest spikeprot.fasta on GISAID
    path = os.path.join(data_root, "spikeprot", "spikeprot.fasta")
    f = open(path, "r")
    f_new = open(path + "_new.fasta", "w")
    count = 0
    for line in f:
        try:
            f_new.write(line.encode('ascii', 'ignore').decode('ascii'))
            count += 1
        except:
            print("Ended")
            break
    
    f.close()
    f_new.close()
    print(count)






## Extract sequences with ID
base_path = os.path.join(data_root, "spikeprot", "spikeprot.fasta_new.fasta")
target_save_path = os.path.join(data_root, "split_data", "{}_record.fasta".format(rbd_name))

def get_ids(filelist):
    Target_ids = []
    for file in filelist:
        with open(os.path.join(csv_path,file), "r") as f:
            for line in f:
                Target_ids.append(line.strip())
    Target_ids = list(set(Target_ids))
    return Target_ids

# read target ids
Target_ids = get_ids(filtered_file)
print("Number of Target:", len(Target_ids))
none_count = 0
target_record = []
total_data = {}
error_record = []

parser = SeqIO.parse(base_path, "fasta")
for j, record in tqdm(enumerate(parser)):
    id = str(record.id)
    try:
        flag = id.split("|")[3]
        total_data[flag] = record
    except:
        error_record.append(record)
for id in tqdm(Target_ids):
    try:
        target_record.append(total_data[id])
        total_data.pop(id)
    except:
        none_count += 1
print("None Number:", none_count)
SeqIO.write(target_record, target_save_path, "fasta")




## Get spike protein and deduplicate
file = os.path.join(data_root, "split_data", "{}_record.fasta".format(rbd_name))
spike = []
error = 0
parser = SeqIO.parse(file, "fasta")
for i, record in tqdm(enumerate(parser)):
    try:
        seq = str(record.seq)
        spike.append(seq)
    except:
        error += 1

spike = list(set(spike))
print(len(spike))
f = open(file + "_spike.txt", "w")
for seq in spike:
    f.write(seq + "\n")
f.close()




## Get RBD sequence
file = os.path.join(data_root, "split_data", "{}_record.fasta_spike.txt".format(rbd_name))
error = 0
counter_rbd=0
f = open(file + "_rbd.txt", "w")
f1 = open(file, "r")
for line in tqdm(f1):
    try:
        seq = line.strip()
        rbd = re.findall("NITN.*KKST", seq)[0]
        if len(rbd) < 256 and len(rbd) > 100:
            f.write(rbd + "\n")
            counter_rbd+=1
    except:
        error += 1
f.close()
f1.close()
print(file, error, counter_rbd)

# deduplicate
del_seq = parent_node[rbd_name]
file = os.path.join(data_root, "split_data", "{}_record.fasta_spike.txt_rbd.txt".format(rbd_name))
count = 0
save_path = os.path.join(data_root, "VOC_train_txt", "rbd_dedup_{}_train.txt".format(rbd_name))
f = open(save_path, "w")
f1 = open(file, "r")
for line in f1:
    if not line.strip() == del_seq:
        f.write(line)
        count += 1
f.close()
f1.close()
print(count)
print("collected sequences saved in {}".format(save_path))




## Calculate frequency
base = parent_node[rbd_name]
total_count = np.zeros(shape=(201,), dtype=np.float64)
base = np.array(list(base))
f = open(save_path, "r")
N = 0
total = 0
for line in tqdm(f):
    total += 1
    qury = line.strip()
    if len(qury) == 201:
        qury = np.array(list(qury))
        count = (base != qury)
        flag_x = np.ones(shape=(201,), dtype=bool)
        flag_x[qury == "X"] = 0
        new = count * flag_x
        if not sum(new) == 0:
            total_count += new
            N += 1
f.close()

probility = total_count / N
new_p = np.zeros(shape=(203,), dtype=np.float64)
new_p[1:-1] = probility
np.save(os.path.join(model_root, "site_mutation_frequency", "{}_probility_mutation_203.npy".format(rbd_name)), new_p)

