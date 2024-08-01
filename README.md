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

[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/Kevinatil/GenPreMut/blob/main/quick_start/SequenceScreen.ipynb)


It will visualize the correctly predicted mutation types in the logo plot.


## Get initial sequence set

To get initial sequence set, please choose one variant as parent node as mutation start variant. For validation experiment, a cutoff node is also needed, which is the endpoint of initial sequence collection.

Please download the latest spike protein sequences on [GISAID](https://www.epicov.org/epi3/frontend), and adjust the data path in `get_init_sequence/get_init_data.py`. If access is denied, it may be necessary to create an account first.

The downloaded fasta file may contain bad lines, which may cause parsing error. If so, please delete bad lines with the following script before parsing.

```python
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
```

The index of the fasta file is formatted as `EPI_ISL_[0-9]+`. The correspondence between fasta index and variant name can be obtained from GISAID.
Run the following script to get the initial sequences, and calculate the mutation frequencies of each site.

```shell
python get_init_sequence/get_init_data.py
```

## PLM finetuning and sequence generation

We adopt [ESM2-650M](https://github.com/facebookresearch/esm) as PLM, and use the collected sequences to adjust the pretrained residue distribution of PLM.

Please adjust the data_path in the script and run as follows.

```shell
python sequence_generate/PLM_finetune/main_ft.py
```

After PLM finetuning, run the following code to generate sequences. The generation scale is set as 1 million. Please adjust the checkpoint path before running.

```shell
python sequence_generate/sequence_generate/generate_sequences.py
```

## Sequence screening

We adopt host-level and herd-level screening. For host-level, we adopt prediction models in [E2VD](https://github.com/ZhiweiNiepku/E2VD) to predict key properties for virus to survive and spread. For herd-level, we propose quantified antibody barrier score to validate the capability of virus to break through the immune barrier of humans.


### Host-level property prediction model

- Data process

We use the dataset from *Deep mutational scans for ACE2 binding, RBD expression, and antibody escape in the SARS-CoV-2 Omicron BA.1 and BA.2 receptor binding domains*.

Run the following script to do sequence filtering and to split training and testing sets.

```shell
python sequence_screening/property_prediction/training/data_process_bind_expr.py
```

- Model training

The training and testing process of binding affinity and expression are the same. Take expression as example.

```shell
python sequence_screening/property_prediction/training/expression_train/expr_train_all.py
```

Calculate the statistical data of training features for test feature normalization.

```shell
python sequence_screening/property_prediction/training/expression_train/statistic.py
```

- Property prediction

Run the predicting script to do expression prediction.

```shell
python sequence_screening/property_prediction/predicting/predict_expr_emb_main.py
```

The predicted results are saved in `.npy` format.

After property prediction, run the merging code to merge prediction results.

```shell
python sequence_screening/ranking/merge_predicts.py
```

### Herd-level quantified antibody barrier score

For the calculation of quantified antibody barrier score, please download the antibody escape data first. We use two antibody escape datasets from two studies.


- [Imprinted SARS-CoV-2 humoral immunity induces convergent Omicron RBD evolution](https://www.nature.com/articles/s41586-022-05644-7)
- [Repeated Omicron infection alleviates SARS-CoV-2 immune imprinting](https://www.nature.com/articles/s41586-023-06753-7)

We first select the antibody data to use and calculate the average antibody escape score of each antibody group.

For antibody barrier model1, please download `use_res_clean.csv` and `supplementary2.csv` at https://github.com/jianfcpku/convergent_RBD_evolution. 

For antibody barrier model2, please download `antibody_dms_merge_no_filter_clean.csv` and `antibody_info.csv` at https://github.com/jianfcpku/SARS-CoV-2-reinfection-DMS.

```shell
python sequence_screening/quantified_antibody_barrier_score/prepare_model1.py

python sequence_screening/quantified_antibody_barrier_score/prepare_model2.py
```

Run the following script to calculate quantified antibody barrier score for each dataset.


### Sequence screening

Run the following script to screen generated sequences, and calculate the frequency of each mutation type.

```shell
python sequence_screening/ranking/screening.py
```

