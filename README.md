# GenPreMut

The official code repository of "Generative prevalent mutation prediction across the entire evolutionary landscape through host-to-herd simulated evolution".

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Repo contents](#repo-contents)
- [Demo and Instructions for Use](#demo-and-instructions-for-use)
- [License](#license)


# Overview

To link selective pressures inside and outside the host, we propose GenPreMut, a generative prediction pipeline for prevalent mutations spanning intra-host evolution to herd-level evolution. The protein language model (PLM) guided by mutation probability are employed as a generator of massive mutations, and the prediction models of selective pressures within and outside the host are adopted as diverse evolution drivers.
A strategy linking intra-host evolution and herd-level evolution is proposed to introduce herd-level selective pressures, transitioning intra-host molecular-level mutational drivers to herd-level antibody barrier.
GenPreMut accurately reproduces previous high-risk mutations, significantly outperforming the state-of-the-art prediction methods.


# System Requirements

## Hardware requirements

The PLM finetuning and downstream predictions using extracted features can be conducted on a standard computer with NVIDIA GPU with no less than 8G memory and enough RAM to run training and inference. The downstream experiments are conducted on Tesla V100 * 1 (32GB).

## Software requirements

### OS Requirements

The experiments are conducted on Ubuntu 16.04.1.

### Python Dependencies

The required python environments are as follows.

```
python>=3.9
torch>=2.0.1
scikit-learn
numpy
seaborn
matplotlib
transformers
datasets
biopython
dmslogo
logomaker
```

Besides, we reproduced the result of [MLAEP](https://github.com/WHan-alter/MLAEP). The required environments can refer to [environment.yml](https://github.com/WHan-alter/MLAEP/blob/master/environment.yml).

# Installation Guide

## Setting up downstream environment

Create a new environment.

```shell
conda create -n gpm python=3.9
conda activate gpm
```

Install pytorch 2.0.1.

```python
pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 --index-url https://download.pytorch.org/whl/cu118
```

Install other dependencies.

```shell
pip install -r requirements.txt
```

It takes less than 10 minutes to install python dependencies on a standard computer.


# Repo contents

```shell
- get_init_sequence/
    - get_init_data.py # collect initial sequence set
- sequence_generate/
    - PLM_finetune/
        - main_ft.py # PLM finetuning
    - sequence_generate/
        - collator.py # data collator for sequence generation
        - generate_sequences.py # sequence generate
- sequence_screening
    - property_prediction/
        - training/ # property prediction model training
        - predicting/ # predict host-level properties
    - quantified_antibody_barrier_score/
        - prepare_model1.py # get cached file for quantified antibody barrier model 1
        - prepare_model2.py # get cached file for quantified antibody barrier model 2
        - antibody_barrier.py # get quantified antibody barrier score
    - ranking/
        - merge_predicts.py # merge host-level prediction results
        - screening.py # sequence screening
- validation_experiments/
    - MLAEP/
    - generate_scale/ # generation scale ablation experiment
        - draw_scale_changing.py
    - rank_change/ # effectiveness validation
        - rank_change.py
    - stability/ # stability validation
        - rank_robust.py
        - screened_num.py
- model_compare/ # compare with other baseline models
    - draw_rank_change.py
    - draw_site_freq.py
    - draw_site_KL.py
- ckpt/ # checkpoints used in experiments
- data/ # data used in experiments
```



# Demo and Instructions for Use

## Quick start

We provide a quick demo for our sequence screening pipeline. BA.2.1 is adopted as parent node.

- PLM finetuning and sequence generation

We provide a sequence subset for quick start, containing 10,000 sequences. This subset is sampled from our collected sequence set for BA.2.1 validation experiment. 

After PLM finetuning, we generate 1,000 sequences based on the adjusted residue distribution.

[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/Kevinatil/GenPreMut/blob/main/quick_start/SequenceGenerate.ipynb)

It takes less than 10 minutes to finetune PLM, and several seconds to generate sequences. Besides, It will visualize the mutation frequencies of the collected sequence set and generated sequences.


- Sequence screening

We provide our predicted properties for 1 million generated sequences with BA.2.1 as parent node. We show the calculation of quantified antibody barrier score, and screen the generated sequences to get the high risk variants as well as the high risk mutation types.



It will visualize the correctly predicted mutation types in the logo plot.


## Get initial sequence set

To get initial sequence set, please choose one variant as parent node as mutation start variant. For validation experiment, a cutoff node is also needed, which is the endpoint of initial sequence collection.

For example, we take BA.2.1 as parent node, and collect sequences after parent node and before BA.3.

<!-- 【爬虫过程怎么写？】 -->

Please download the latest spike protein sequences on [GISAID](https://www.epicov.org/epi3/frontend), and adjust the data path in `get_init_sequence/get_init_data.py`. If access is denied, it may be necessary to create an account first.

The downloaded fasta file may contain bad lines, which may cause parsing error. If so, please delete bad lines with the following script before parsing.

```python
path = os.path.join(data_root, "spikeprot", "spikeprot0304.fasta")
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
```

Run the following script to get the initial sequences, and calculate the mutation frequencies of each site.

```shell
python get_init_sequence/get_init_data.py
```

