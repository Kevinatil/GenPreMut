import os
import re
import matplotlib.pyplot as plt

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer, AutoModelForMaskedLM

from collator import DataCollatorForMaskedGeneration

rbd_dict={
    'BA.2.1':   'NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST',
    'BA.5.1':   'NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVGGNYNYRYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGVNCYFPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST',
    'XBB.1.5':  'NITNLCPFHEVFNATTFASVYAWNRKRISNCVADYSVIYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKPSGNYNYLYRLFRKSKLKPFERDISTEIYQAGNKPCNGVAGPNCYSPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST',
}
mutation_dict={
    'BA.2.1':   ['G339D','S371F','S373P','S375F','T376A','D405N','R408S','K417N','N440K','S477N','T478K','E484A','Q493R','Q498R','N501Y','Y505H'],
    'BA.5.1':   ['G339D','S371F','S373P','S375F','T376A','D405N','R408S','K417N','N440K','L452R','S477N','T478K','E484A','F486V','Q498R','N501Y','Y505H'],
    'XBB.1.5':  ['G339H','R346T','L368I','S371F','S373P','S375F','T376A','D405N','R408S','K417N','N440K','V445P','G446S','N460K','S477N','T478K','E484A','F486P','F490S','Q498R','N501Y','Y505H'],
}

def get_index(line):
    return int(re.findall(r'[A-Z]([0-9]+)[A-Z]',line)[0])
def get_mut(line):
    return re.findall(r'[A-Z][0-9]+([A-Z])',line)[0]

def numpy_mask_tokens(inputs, probility_mutation, mask_token_id):
    """
    Prepare masked tokens inputs/labels for masked language modeling: 80% MASK, 10% random, 10% original.
    """
    # Numpy doesn't have bernoulli, so we use a binomial with 1 trial
    masked_indices = np.random.binomial(1, probility_mutation, size=probility_mutation.shape).astype(bool)
    masked_lm_positions = np.where(masked_indices == True)[0]
    # The rest of the time (10% of the time) we keep the masked input tokens unchanged
    inputs[masked_lm_positions] = mask_token_id
    return inputs, masked_lm_positions


class GenerateDataset(Dataset):
    def __init__(self, num_samples, rbd_seq, tokenizer):
        super().__init__()
        self.num_samples = num_samples
        self.rbd_seq = rbd_seq
        self.tokenizer = tokenizer

    def __getitem__(self, index):
        return self.tokenizer(self.rbd_seq)

    def __len__(self):
        return self.num_samples

def dump(seqs, name, path):
    file_name = os.path.join(path, "mutation_{}_{}.txt".format(rbd_name, name))
    print('seq num: {}, save path: {}'.format(len(seqs), file_name))
    f = open(file_name, "w")
    for seq in seqs:
        f.write(seq + "\n")
    f.close()


def _check(path):
    raw = rbd_dict['BA.2.1']
    nums = np.zeros(len(raw))
    seqs = set()
    for i in range(10):
        path_ = os.path.join(path, 'mutation_BA.2.1_{}.txt'.format(i))
        f = open(path_, 'r')
        for line in f:
            nums += (np.array(list(line.strip())) != np.array(list(raw)))
            seqs.add(line.strip())
    print(nums)
    print(len(seqs))
    plt.plot(nums)
    plt.savefig('plot.png')


if __name__ == "__main__":
    device = torch.device("cuda")

    rbd_name = 'BA.2.1' # 'BA.5.1' # 'XBB.1.5'
    ft_path = 'ckpt/finetune/checkpoint-xx'
    data_root = 'data/'
    site_freq_path = 'ckpt/site_mutation_frequency/{}_mutation_frequency_203.npy'.format(rbd_name)
    mutation_save_path = os.path.join(data_root, 'raw_seqs')

    total_number = 1_000_000
    step = 10_000

    tokenizer = AutoTokenizer.from_pretrained(ft_path)
    model = AutoModelForMaskedLM.from_pretrained(ft_path).to(device).eval()

    max_len = 203
    max_mask = 5
    topk = 10
    batch_size = 4

    rbd_seq = rbd_dict[rbd_name]
    rbd_id = np.array(tokenizer(rbd_seq)['input_ids'])

    save_steps=np.arange(0, total_number+1, step)[1:]
    os.makedirs(mutation_save_path, exist_ok=True)

    probility_mutation = np.load(site_freq_path)

    collator = DataCollatorForMaskedGeneration(tokenizer, torch.tensor(probility_mutation), max_mask, device=device)

    dataset = GenerateDataset(num_samples=total_number*100, 
                              rbd_seq=rbd_seq, 
                              tokenizer=tokenizer)
    dataloader = DataLoader(dataset, batch_size=batch_size, collate_fn=collator)

    with torch.no_grad():
        output_seq = set()
        output_seq_tp = set()
        process = 0
        for i, data in enumerate(dataloader):
            if i % 500 == 0:
                print('>>>>> {} loops, current sequence num: {}'.format(i, len(output_seq)))
            token_ids_ = data['token_ids'].cpu()
            data.pop('token_ids')

            out = model(**data)
            indices = torch.topk(out['logits'], topk, dim = -1).indices.cpu()
            bs = indices.shape[0]
            indices = indices.reshape(-1, topk)
            for _ in range(20):
                token_ids = token_ids_.clone()
                index_ran = np.random.randint(0, topk, size=(indices.shape[0]))
                predict_id = indices[range(indices.shape[0]), index_ran]
                predict_id = predict_id.reshape(bs, -1)

                mask_ = (token_ids == tokenizer.mask_token_id)&(predict_id >= 4)&(predict_id <= 23)
                token_ids[mask_] = predict_id[mask_]

                sequences = tokenizer.batch_decode(token_ids, skip_special_tokens=True, clean_up_tokenization_spaces=True)

                for sequence in sequences:
                    sequence = re.sub(r'\s', '', sequence)
                    if sequence not in output_seq:
                        output_seq.add(sequence)
                        output_seq_tp.add(sequence)

                    if process >= len(save_steps):
                        break
                    if len(output_seq) >= save_steps[process]:
                        print('process {}, output_seq: {}, output_seq_tp: {}'.format(process, len(output_seq), len(output_seq_tp)))
                        dump(output_seq_tp, process, mutation_save_path)
                        process+=1
                        output_seq_tp=set()

                if process >= len(save_steps):
                    break
            if process >= len(save_steps):
                break
        print(total_number, len(output_seq))
