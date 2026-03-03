import argparse
from utils import seed_everything, get_device, FocalLoss
from dataset import COData
from trainer import Trainer
import os
import torch
import json
import torch.nn as nn   
from dicts_64 import genetic_dictionary_64
from dicts_12 import genetic_dictionary_12
from dicts_7 import genetic_dictionary_7

from Bio import SeqIO
import pandas as pd


def convert_to_python_types(obj):
    """
    Recursively convert tensors in a nested structure (like dict or list) to standard Python types.
    """
    if isinstance(obj, dict):
        return {k: convert_to_python_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_python_types(v) for v in obj]
    elif isinstance(obj, torch.Tensor):
        return obj.item() if obj.numel() == 1 else obj.tolist()
    else:
        return obj

def run(args_dict: dict) -> None:
    """
    Run the training process based on the provided arguments.
    """
    seed_everything(args_dict['seed'])
    device = get_device(args_dict['device'])
    
    if args_dict['dict_version'] == 64:
        genetic_dictionary = genetic_dictionary_64
    elif args_dict['dict_version'] == 12:
        genetic_dictionary = genetic_dictionary_12
    elif args_dict['dict_version'] == 7:
        genetic_dictionary = genetic_dictionary_7
    else:
        raise ValueError("Dictionary version not recognized.")

    loss_fn = nn.CrossEntropyLoss() if args_dict['loss_fn'] == 'ce' else FocalLoss()

    file_path = str(args_dict['data_name']) + '/' + str(args_dict['species']) + '/' + str(args_dict['name']) + '/' + str(args_dict['file_name']) + '/' + str(args_dict['seed']) + '/'
    trainer = Trainer(
        seed=args_dict['seed'],
        device=device,
        file_name=file_path,
        loss_fn=loss_fn,
        genetic_dictionary=genetic_dictionary
    )

    data = COData(
        data_name=args_dict['data_name'],
        species=args_dict['species'], 
        clip_size=args_dict['clip_size'], 
        seed=args_dict['seed'],
        genetic_dictionary=genetic_dictionary,
    )
    data.load_data()
    trainer._set_data(data, hgcOPT=args_dict['hgcOPT'])
    trainer.define_model(args_dict)

    if args_dict['epoch_pt'] > 0:
        trainer.set_loaders_pt() # size_info=[[0, 100000, 512, 256],])
        trainer.define_pretraining_modules()
        trainer.pretraining(
            epochs=args_dict['epoch_pt'],
            optimizer=args_dict['optimizer_pt'],
            lr=args_dict['lr_pt'],
            verbose=args_dict['verbose'],
            mlm_codon=args_dict['mlm_codon'],
            mlm_aa=args_dict['mlm_aa'],
            contrastive_pt=args_dict['contrastive_pt'],
        )

    trainer.set_loaders() # size_info=[[0, 100000, 256, 512]]

    trainer.fit(
        epochs=args_dict['epoch'],
        optimizer=args_dict['optimizer'],
        lr=args_dict['lr'],
        earlystop=args_dict['earlystop'],
        wo_pad=args_dict['wo_pad'],
        verbose=args_dict['verbose'],
        lambda_codon=args_dict['lambda_codon'],
        lambda_emb=args_dict['lambda_emb'],
        lambda_mfe=args_dict['lambda_mfe'],
        mcdropout_rate=args_dict['mcdropout_rate'],
        iter_mcdropout=args_dict['iter_mcdropout'],
    )
    converted_logs_fit = convert_to_python_types(trainer.log_fit)

    log_path='./logs/'+ file_path
    os.makedirs(log_path, exist_ok=True)
    with open(log_path + 'logs.json', 'w') as file:
        json.dump(converted_logs_fit, file, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Train a model based on provided arguments.")

    parser.add_argument('--data_name', type=str, default='all', choices=['all', 'genart', 'gensmart', 'novopro', '3000yc', 'linear3000'])
    parser.add_argument('--species', type=str, default='homosapiens', choices=['ecoli', 'homosapiens', 'yeast', 'musculus'])
    parser.add_argument('--clip_size', type=int, default=512, help='Clip size for sequences.')
    parser.add_argument('--dict_version', type=int, default=7, choices=[64, 12, 7], help='Dictionary version.')
    parser.add_argument('--hgcOPT', action='store_true', help='Use hgcOPT dataset.')
    
    # Model arguments
    parser.add_argument('--name', type=str, default='coformer', choices=['coformer', 'lstmsimple', 'realjmformer', 'lstm6', 'lstm4', 'transformer' ], help='Model name.')
    parser.add_argument('--depth', type=int, default=3, help='Depth (layers) of backbone model.')
    parser.add_argument('--dim_in', type=int, default=32, help='input embedding dimension')
    parser.add_argument('--dim_attn', type=int, default=32, help='attention embedding dimension')
    parser.add_argument('--dim_out', type=int, default=12, help='output embedding dimension')

    parser.add_argument('--n_heads', type=int, default=4, help='number of heads')
    parser.add_argument('--dropout_attn', type=int, default=0.1, help='dropout rate for attention')

    parser.add_argument('--mult_ff', type=int, default=2, help='multiplication for feedforward')
    parser.add_argument('--dropout_ff', type=int, default=0.1, help='dropout rate for feedforward')

    # Training arguments
    parser.add_argument('--epoch', type=int, default=200, help='Number of training epochs.')
    parser.add_argument('--batchsize', type=int, default=512, help='Training batch size.')
    parser.add_argument('--optimizer', type=str, default='adamw', choices=['adamw', 'adam', 'sgd'], help='Optimizer choice.')
    parser.add_argument('--lr', type=float, default=1e-4, help='Learning rate.')
    parser.add_argument('--earlystop', type=int, default=30, help='Early stopping.')
    parser.add_argument('--wo_pad', action='store_true')
    parser.add_argument('--verbose', action='store_false')
    parser.add_argument('--file_name', type=str, default='best_model', help='file name for saving the best model.')
    parser.add_argument('--loss_fn', type=str, default='ce', choices=['ce', 'focal'], help='Loss function choice.')

    parser.add_argument('--lambda_codon', type=float, default=0.0, help='Lambda for codon loss.')
    parser.add_argument('--lambda_emb', type=float, default=0.0, help='Lambda for embedding loss.')
    parser.add_argument('--lambda_mfe', type=float, default=0.0, help='Lambda for mfe loss.')
    parser.add_argument('--window', type=int, default=64, help='Window size for ce loss.')
    parser.add_argument('--step_size', type=int, default=32, help='Step size for ce loss.')


    parser.add_argument('--iter_mcdropout', type=int, default=0, help='Number of samples for Monte Carlo dropout.')
    parser.add_argument('--mcdropout_rate', type=float, default=0.0, help='Dropout rate for Monte Carlo dropout.')

    # Pretrainint arguments
    parser.add_argument('--epoch_pt', type=int, default=0, help='Number of training epochs.')
    parser.add_argument('--batchsize_pt', type=int, default=128, help='Training batch size.')
    parser.add_argument('--optimizer_pt', type=str, default='adamw', choices=['adamw', 'adam', 'sgd'], help='Optimizer choice.')
    parser.add_argument('--lr_pt', type=float, default=1e-3, help='Learning rate.')
    parser.add_argument('--mlm_codon', action='store_true')
    parser.add_argument('--mlm_aa', action='store_true')
    parser.add_argument('--contrastive_pt', action='store_true')
    
    # Misc arguments
    parser.add_argument('--seed', type=int, default=1, help='Random seed for reproducibility.')
    parser.add_argument('--device', type=int, default=-1, choices=[-1, 0, 1], help='Device choice (-1 for CPU, 0 or 1 for GPU).')
    parser.add_argument('--log_name', type=str, default='defalut', help='log name for saving the logs.')

    args = parser.parse_args()
    args_dict = vars(args)
    run(args_dict)