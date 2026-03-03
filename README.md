# Codon Optimization Framework

A deep learning framework for multi-objective codon optimization integrating codon-level accuracy, embedding consistency, and mRNA stability (MFE) regularization.

---

## 🧬 Overview

This project provides a flexible training pipeline for codon optimization across multiple species, supporting:

- Multi-objective learning (Codon + Embedding + MFE)
- Transformer/LSTM-based architectures
- Masked language modeling (MLM) pretraining
- Contrastive pretraining
- Monte Carlo dropout-based uncertainty estimation

Supported species:
- `homosapiens`
- `ecoli`
- `yeast`
- `musculus`

---

## 📂 Project Structure

```
project/
 ├── data/
 │    └── homosapiens/
 ├── models/
 ├── src/
 ├── logs/
 ├── train.py
 └── README.md
```

---

## ⚙ Installation

```bash
git clone https://github.com/username/Codon_Optimization.git
cd Codon_Optimization
pip install -r requirements.txt
```

---

## 🚀 Quick Start

Basic training example:

```bash
python train.py \
    --data_name all \
    --species homosapiens \
    --name coformer \
    --epoch 200 \
    --batchsize 512 \
    --lambda_codon 1.0 \
    --lambda_mfe 0.5
```

---

# 🧠 Argument Description

---

## 🔹 Dataset Arguments

| Argument | Description |
|----------|-------------|
| `--data_name` | Dataset source (`all`, `genart`, `gensmart`, `novopro`, `3000yc`, `linear3000`) |
| `--species` | Target species (`ecoli`, `homosapiens`, `yeast`, `musculus`) |
| `--clip_size` | Maximum sequence length |
| `--dict_version` | Codon dictionary version (`64`, `12`, `7`) |
| `--hgcOPT` | Use HGC optimization dataset |

---

## 🔹 Model Arguments

| Argument | Description |
|----------|-------------|
| `--name` | Model architecture (`coformer`, `lstmsimple`, `realjmformer`, `lstm6`, `lstm4`, `transformer`) |
| `--depth` | Number of backbone layers |
| `--dim_in` | Input embedding dimension |
| `--dim_attn` | Attention embedding dimension |
| `--dim_out` | Output embedding dimension |
| `--n_heads` | Number of attention heads |
| `--dropout_attn` | Attention dropout rate |
| `--mult_ff` | Feed-forward expansion factor |
| `--dropout_ff` | Feed-forward dropout rate |

---

## 🔹 Training Arguments

| Argument | Description |
|----------|-------------|
| `--epoch` | Training epochs |
| `--batchsize` | Batch size |
| `--optimizer` | Optimizer (`adamw`, `adam`, `sgd`) |
| `--lr` | Learning rate |
| `--earlystop` | Early stopping patience |
| `--loss_fn` | Loss function (`ce`, `focal`) |
| `--lambda_codon` | Codon loss weight |
| `--lambda_emb` | Embedding loss weight |
| `--lambda_mfe` | MFE loss weight |
| `--window` | Sliding window size |
| `--step_size` | Step size for window |
| `--iter_mcdropout` | Number of MC dropout samples |
| `--mcdropout_rate` | MC dropout rate |

---

## 🔹 Pretraining Arguments

| Argument | Description |
|----------|-------------|
| `--epoch_pt` | Pretraining epochs |
| `--batchsize_pt` | Pretraining batch size |
| `--optimizer_pt` | Pretraining optimizer |
| `--lr_pt` | Pretraining learning rate |
| `--mlm_codon` | Enable codon-level MLM |
| `--mlm_aa` | Enable amino acid MLM |
| `--contrastive_pt` | Enable contrastive pretraining |

---

## 🔹 Miscellaneous

| Argument | Description |
|----------|-------------|
| `--seed` | Random seed |
| `--device` | Device (`-1` CPU, `0` or `1` GPU) |
| `--file_name` | Best model filename |
| `--log_name` | Log filename |

---

## 📊 Objective Function

The total loss can be expressed as:

```
L_total = L_ce
        + λ_codon * L_codon
        + λ_emb * L_embedding
        + λ_mfe * L_mfe
```

---

## 🧪 Example Configurations

**Codon-only training**

```bash
python train.py --lambda_codon 1.0 --lambda_mfe 0.0
```

**Multi-objective (Codon + MFE)**

```bash
python train.py --lambda_codon 1.0 --lambda_mfe 0.5
```

**With contrastive pretraining**

```bash
python train.py --contrastive_pt --epoch_pt 50
```

---

## 📌 Reproducibility

```bash
python train.py --seed 1
```

---

## 📜 License

Specify your license here.

---

## 📖 Citation

If you use this work, please cite:

```
Your paper citation here
```
