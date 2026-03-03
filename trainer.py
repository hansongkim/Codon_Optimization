# from math import inf
from torch.utils.data import DataLoader
from models.modelwrapper import ModelWrapper, Augmentation
from utils import get_optimizer, FocalLoss
from dataset import COData, CODataset
from optimizer import Optimizer
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, precision_score, recall_score, matthews_corrcoef
from dicts_64 import index2codon_wopad, index2aminoacid 
from dicts_7 import index2codon_global, aa2global_index, codon_id2triplet, codon_to_aa
from torch.nn.utils.rnn import pad_sequence
from sklearn.cluster import KMeans

import copy
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import pandas as pd
import os
import json
from Bio import SeqIO

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from sklearn.manifold import TSNE
import umap.umap_ as umap
import os
import glob



from collections import Counter


class Trainer:
    def __init__(
            self, 
            species='homosapiens',
            genetic_dictionary=None,
            seed=1, 
            device='cpu', 
            file_name='best_model', 
            loss_fn=nn.CrossEntropyLoss()
    ):
        self.seed = seed
        self.device = device
        self._setup_seed()
        self.loss_fn = loss_fn
        self.optim = None
        self.best_model_state = None
        if genetic_dictionary is None:
            raise ValueError("Genetic dictionary must be provided.")
        self.optimizer = Optimizer(seed=seed, device=device, species=species, genetic_dictionary=genetic_dictionary)
        self._set_file_paths(file_name)
        self.verbose = True
        self.ticks=['\\', '-', '/', '|']

    def _setup_seed(self):
        """Sets the seed for reproducibility."""
        torch.manual_seed(self.seed)
        np.random.seed(self.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(self.seed)
    
    def _create_directory(self, path):
        """Create directory if it does not exist."""
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"Created directory: {path}")
        else:
            print(f"Directory already exists: {path}")

    def _set_file_paths(self, file_name):
        """Set file paths for saving and loading models."""
        self.dir_path = f'./saved_models/{file_name}/'
        self._create_directory(self.dir_path)
        
        self.log_path = f'./logs/{file_name}/'
        self._create_directory(self.log_path)
        
        # Model
        self.model_state_path = f'{self.dir_path}model_state.pth'
        self.model_path = f'{self.dir_path}model.pth'
        # Embedding
        self.embedding_state_path = f'{self.dir_path}embedding_state.pth'
        self.embedding_path = f'{self.dir_path}embedding.pth'
        # Backbone
        self.backbone_state_path = f'{self.dir_path}backbone_state.pth'
        self.backbone_path = f'{self.dir_path}backbone_state.pth'
        # To out
        self.to_out_state_path = f'{self.dir_path}to_out_state.pth'
        self.to_out_path = f'{self.dir_path}to_out.pth'
        # pretrained
        self.model_state_pt_path = f'{self.dir_path}model_pt_state.pth'
        self.model_pt_path = f'{self.dir_path}model_pt.pth'
        # Parameters
        self.model_params_path = f'{self.dir_path}params.json'

    def _save_params(self, params): # Save the model parameters
        with open(self.model_params_path, 'w') as file:
            json.dump(params, file, indent=4)

    def split_dataset(self, df, split_rate=[0.6, 0.2, 0.2]):
        list_entry = df['Entry'].unique()
        list_train, tmp_data = train_test_split(list_entry, test_size=split_rate[1]+split_rate[2], random_state=self.seed)    
        list_valid, list_test = train_test_split(tmp_data, test_size=0.5, random_state=self.seed)
        df_train = df[df['Entry'].isin(list_train)]
        df_valid = df[df['Entry'].isin(list_valid)]
        df_test = df[df['Entry'].isin(list_test)]
        return df_train, df_valid, df_test
        
    def _set_data(self, data, split_rate=[0.6, 0.2, 0.2], hgcOPT=False):
        self.data = data
        self.data.split_dataset(split_rate)
        self.df_train, self.df_valid, self.df_test = self.data.df_train, self.data.df_valid, self.data.df_test

        self.df_test.to_csv("test_set.csv", index=False)
        
        if hgcOPT:
            file_path = "./hgcOPT.fa"
            for sample in SeqIO.parse(file_path, "fasta"):
                name = sample.name
                dna = sample.seq
                rna = sample.seq.transcribe()
                protein = str(rna.translate())
                aminoacid = [a for a in protein]
                codon = [str(rna[i:i+3]) for i in range(0, len(rna), 3)]
            seq_dna = ''
            for i in dna:
                seq_dna+=i

            seq_rna = ''
            for i in rna:
                seq_rna+=i   
            codon = data.seq2codon(seq_rna)
            self.data.df_hgcOPT = pd.DataFrame([[name, seq_dna, seq_rna, aminoacid, len(aminoacid), codon, len(codon), len(codon)]],  
                                               columns=['Name', 'DNA', 'RNA', 'AminoAcid', 'length', 'Codon', 'length_codon', 'length_codon_wopad'])
        else:
            self.data.df_hgcOPT = None

    def set_loaders(self, size_info=None):
        """Sets up the DataLoader for training, validation, and testing."""
        if size_info is None: # list for pairs of ['start', 'end', 'clip_size', 'batch_size']
            size_info = [ [0, 128, 128, 256],
                         [128, 256, 256, 128],
                         [256, 512, 512, 64],
                         [512, 100000, 1024, 32]]
        self.train_loader =[]
        self.valid_loader = []
        self.test_loader = []

        for info in size_info:
            try:
                print(f"Setting up DataLoader for length: {info[0]} - {info[1]}")
                df_train_tmp = self.df_train[(self.df_train['length'] > info[0]) & (self.df_train['length'] <= info[1])]
                df_valid_tmp = self.df_valid[(self.df_valid['length'] > info[0]) & (self.df_valid['length'] <= info[1])]
                df_test_tmp = self.df_test[(self.df_test['length'] > info[0]) & (self.df_test['length'] <= info[1])]

                data_train, data_valid, data_test = self.data.prepare_train_valid_test(df_train_tmp, df_valid_tmp, df_test_tmp, clip_size=info[2])

                # data : [aa, codon, label, pads] 
                # CODataset : [aa, input(aa or codon), label, pad]
                self.train_loader.append(DataLoader(CODataset(data_train[0], data_train[1], data_train[2], data_train[3]), info[3], shuffle=True))
                self.valid_loader.append(DataLoader(CODataset(data_valid[0], data_valid[1], data_valid[2], data_valid[3]), info[3], shuffle=False))
                self.test_loader.append(DataLoader(CODataset(data_test[0], data_test[1], data_test[2], data_test[3]), info[3], shuffle=False))
            except:
                print(f"Data for length {info[0]} - {info[1]} not found.")
                continue

    def define_model(self, params, len_aminoacids=None, len_codons=None, p_mask=0.4):
        if len_aminoacids is None or len_codons is None:
            params['len_aminoacids'] = self.data.len_aminoacids
            params['len_codons'] = self.data.len_codons
            self.len_aminoacids = params['len_aminoacids']
            self.len_codons = params['len_codons']

        params['dim_out'] = self.data.len_labels
        self.dim_out = params['dim_out']
        self.model = ModelWrapper(params).to(self.device)
        self.best_model_state = self.model.state_dict()
        self.augmentation = Augmentation(dim=self.model.params['dim_in'], p_mask=p_mask).to(self.device)
        self._save_params(params)

    def define_pretraining_modules(self, p_mask=0.4):
        if self.model is None:
            raise ValueError("Model must be defined first. Use define_model() method.")
        self.model._set_pretraining_modules(device=self.device)
        self.best_model_pt_state = self.model.state_dict()
    
    def reset_pretraining_modeuls(self):
        self.model.to_codon = None
        self.model.to_aa = None
        self.model.to_contrast = None
        self.augmentation = None
        
    def save_model(self):
        torch.save(self.model, self.model_path)
        torch.save(self.model.embedding, self.embedding_path)
        torch.save(self.model.backbone, self.backbone_path)
        torch.save(self.model.to_out, self.to_out_path)

    def load_model(self):
        self.model = torch.load(self.model_path)
        self.model.eval()
        self.model.embedding = torch.load(self.embedding_path)
        self.model.embedding.eval()
        self.model.backbone = torch.load(self.backbone_path)
        self.model.backbone.eval()
        self.model.to_out = torch.load(self.to_out_path)
        self.model.to_out.eval()

    def save_module(self, module, path):
        torch.save(module, path)

    def load_module(self, path):
        self.model = torch.load(path)
        self.model.eval()

    def save_model_state(self, path):
        self.best_model_state = copy.deepcopy(self.model.state_dict())
        self.best_embedding_state = copy.deepcopy(self.model.embedding.state_dict())
        self.best_backbone_state = copy.deepcopy(self.model.backbone.state_dict())
        self.best_to_out_state = copy.deepcopy(self.model.to_out.state_dict())

        torch.save(self.best_model_state, path)
        torch.save(self.best_embedding_state, self.embedding_state_path)
        torch.save(self.best_backbone_state, self.backbone_state_path)
        torch.save(self.best_to_out_state, self.to_out_state_path)

    def load_model_state(self, module_name=None, path=None):
        if module_name == 'embedding':
            self.model.embedding.load_state_dict(torch.load(path), strict=False)
        elif module_name == 'backbone':
            self.model.backbone.load_state_dict(torch.load(path), strict=False)
        elif module_name == 'to_out':
            self.model.to_out.load_state_dict(torch.load(path), strict=False)
        else:
            if path is None:
                path = self.model_state_path
            self.model.load_state_dict(torch.load(path), strict=False)

    def set_mcdropout(self, n_iter=0, rate=0.0):
        self.mcdropout = nn.Dropout(p=0.0)
        self.optimizer.mcdropout = nn.Dropout(p=0.0)
        self.iter_mcdropout = n_iter
        self.optimizer.iter_mcdropout = self.iter_mcdropout
        if n_iter > 0:
            self.mcdropout = nn.Dropout(p=rate)
            self.optimizer.mcdropout = nn.Dropout(p=rate)
   
    def logging_loss(self, losses, verbose=True):
        if verbose:
            print('\r|  total  |aminoacid|  codon  |  cosine |   mfe   |')
            print(f'\r| {losses[0]:.5f} | {losses[1]:.5f} | {losses[2]:.5f} | {losses[3]:.5f} | {losses[4]:.5f} |')
        return {'total': losses[0], 'aa': losses[1], 'codon': losses[2], 'cosine': losses[3], 'mfe': losses[4]}
    
    def logging_score(self, scores, stats, verbose=True):
        if verbose:
            print('\r| accuracy| f1 macro|f1 weight|   mcc   |   CAI   |  GC   |   GC3   |')
            print(f'\r| {scores[0]:.5f} | {scores[1]:.5f} | {scores[2]:.5f} | {scores[3]:.5f} | {stats[0]:.5f} | {stats[1]:.5f} | {stats[2]:.5f} |')
        return {'acc': scores[0], 'f1_macro': scores[1], 'f1_weight': scores[2], 'mcc': scores[3], 'CAI': stats[0], 'GC': stats[1],  'GC3': stats[2]}

    def fit(
            self, 
            epochs=100, 
            optimizer='adamw', 
            lr=1e-3, 
            earlystop=15, 
            wo_pad=False, 
            verbose=True, 
            lambda_codon=0.0,
            lambda_emb=0.0,
            lambda_mfe=0.0,
            window=16,
            step_size=8,
            mcdropout_rate=0.0,
            iter_mcdropout=0
        ):
        self.log_fit = {'train': {}, 'valid': {}, 'test': {}}
        self.epochs = epochs
        self.wo_pad = wo_pad
        self.verboes = verbose
        self.lambda_codon = lambda_codon
        self.lambda_emb = lambda_emb
        self.lambda_mfe = lambda_mfe
        self.window = window
        self.step_size = step_size
        
        islast = False

        self.set_mcdropout(iter_mcdropout, mcdropout_rate)
        self.optim = get_optimizer(params=self.model.parameters(), opt_name=optimizer, lr=lr)

        earlystop = earlystop if earlystop is not None else epochs
        best_valid, patience = 0, 0


        self.train_mode = False
        
        if self.data.df_hgcOPT is not None:
            stat_hgcOPT = self.get_hgcOPT_statistics(self.data.df_hgcOPT)
            if self.verbose:
                print('| Data  |Epoch|   CAI   |  GC   |   GC3   |   MFE  |')
                print(f'| hgcOPT |  {0}  | {stat_hgcOPT[0]:.5f} | {stat_hgcOPT[1]:.5f} | {stat_hgcOPT[2]:.5f} | {stat_hgcOPT[3]:.5f} |')
                print(F'Optimized RNA sequence: {stat_hgcOPT[4]}')
                print(F'Optimized codon sequence: {stat_hgcOPT[5]}')

        for epoch in range(1, epochs + 1):

            islast = False
            
            self.train_mode = True
            self.model.train()
            l_train, s_train = self.step_epoch(self.train_loader, islast)

            if self.verbose:
                print(f'\r[Train][{epoch}/{epochs}]_______________________________________________________________')
            self.log_fit['train'][epoch] = {'losses': self.logging_loss(l_train)}

            self.train_mode = False
            self.model.eval()
            with torch.no_grad():
                l_valid, s_valid = self.step_epoch(self.valid_loader, islast)
                stat_valid= self.get_statistics(self.data.df_valid)
                if self.verbose:
                    print(f'\r[Valid][{epoch}/{epochs}]______________________________________________________________')
                self.log_fit['valid'][epoch] = {'losses': self.logging_loss(l_valid, verbose=False), 'scores': self.logging_score(s_valid, stat_valid)}
                

                score_valid = s_valid[0]
                isBest = True if score_valid > best_valid else False 
                if isBest:
                    
                    islast = True
                    
                    e_best, patience = epoch, 0
                    best_valid = score_valid
                    l_b_valid, s_b_valid = l_valid, s_valid
                    l_test, s_test = self.step_epoch(self.test_loader, islast)
                    stat_test= self.get_statistics(self.data.df_test)





                    #                     # =====================================
                    # # ✅ Embedding visualization (best model)
                    # # =====================================
                    # with torch.no_grad():
                    #     all_ctx, all_seq = [], []
                    #     for batch in self.test_loader:
                    #         seq_aa = batch[0].to(self.device)
                    #         seq_label = batch[2].to(self.device)
                    #         pad = batch[3].to(self.device)
                            
                    #         seq_pred, context_aa = self.forward(seq_aa)
                    #         all_ctx.append(context_aa.detach().cpu())
                    #         all_seq.append(seq_aa.detach().cpu())
                
                    #     all_ctx = torch.cat(all_ctx, dim=0)
                    #     all_seq = torch.cat(all_seq, dim=0)
                
                    #     # 실제 시각화 호출 (로컬 매핑 7~27 버전)
                    #     self.visualize_embedding(all_ctx, all_seq)


                    # if islast:
                    #     self.last_context_aa = torch.cat(all_context, dim=0)
                    #     self.last_seq_aa = torch.cat(all_seq, dim=0)
                    #     self.step_epoch(self.test_loader, islast)
                    #     self.visualize_embedding(self.last_context_aa, self.last_seq_aa)






                    

                    if self.data.df_hgcOPT is not None:
                        stat_hgcOPT = self.get_hgcOPT_statistics(self.data.df_hgcOPT)
                    if self.verbose:
                        print(f'\r[test][{epoch}/{epochs}]_______________________________________________________________')
                    self.log_fit['test'][e_best] = {'losses': self.logging_loss(l_test, verbose=False), 'scores': self.logging_score(s_test, stat_test)}
                    self.save_model_state(self.model_state_path)
                patience += 1
                                                     
            if verbose:
                print(f'---------------------------------------- [Best] ----------------------------------------')
                print('| Data  |Epoch| | loss | accuracy | f1 macro|f1 weight|   mcc   |   CAI   |  GC   |   GC3   |   MFE  |')
                print(f'| Valid |  {e_best}  | {l_b_valid[0]:.5f} | {s_b_valid[0]:.5f} | {s_b_valid[1]:.5f} | {s_b_valid[2]:.5f} |', end='')
                print(f' {s_b_valid[3]:.5f} | {stat_valid[0]:.5f} | {stat_valid[1]:.5f} | {stat_valid[2]:.5f} | {0} |')
                print(f'| test  |  {e_best}  | {l_test[0]:.5f} | {s_test[0]:.5f} | {s_test[1]:.5f} | {s_test[2]:.5f} |', end='')
                print(f' {s_test[3]:.5f} | {stat_test[0]:.5f} | {stat_test[1]:.5f} | {stat_test[2]:.5f} | {0} |')
                if self.data.df_hgcOPT is not None:
                    print('| Data  |Epoch|   CAI   |  GC   |   GC3   |   MFE  |')
                    print(f'| hgcOPT |  {e_best}  | {stat_hgcOPT[0]:.5f} | {stat_hgcOPT[1]:.5f} | {stat_hgcOPT[2]:.5f} | {stat_hgcOPT[3]:.5f} |')
                    print(F'Optimized RNA sequence: {stat_hgcOPT[4]}')
                    print(F'Optimized codon sequence: {stat_hgcOPT[5]}')
                print(f'---------------------------------------------------------------------------------------')
                print()


            # # 여기 수정하기 ################################
            # if epoch % 1 == 0:
            #     with torch.no_grad():
            #         len_aa = self.model.params['len_aminoacids']
            #         aa_labels = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            #                      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
            #         embedding_tensor = self.model.embedding.weight[:len_aa].detach().cpu().numpy()
            
            #         reducer = umap.UMAP(n_neighbors=10, min_dist=0.1, random_state=42)
            #         reduced = reducer.fit_transform(embedding_tensor)
            
            #         os.makedirs('embedding_plots', exist_ok=True)
            #         plt.figure(figsize=(8, 6))
            #         for i, label in enumerate(aa_labels):
            #             x, y = reduced[i]
            #             plt.scatter(x, y, s=50)
            #             plt.text(x + 0.2, y + 0.2, label, fontsize=12)
            #         plt.title(f'Amino Acid Embedding - Epoch {epoch}')
            #         plt.grid(True)
            #         plt.axis('off')
            #         plt.tight_layout()
            #         plt.savefig(f'embedding_plots/umap_embedding_epoch_{epoch}.png', dpi=300)
            #         plt.close()
            #         print(f"[✅] Embedding plot saved for epoch {epoch}")
            #     # print("aa embedding?:", embedding_tensor)
            #     # print("embedding shape?: ", embedding_tensor.shape)

            # print("index2aminoacid in trainer:", index2aminoacid)
            # if patience > earlystop:
            #     if verbose:
            #         print('Early-stopping....!')
            #     break

        self.load_model_state()



    def visualize_embedding_codon_vocab(self, context_codon, seq_codon):
        """
        seq_codon (vocab ID: 22~85)을 기준으로
        codon 임베딩을 t-SNE로 시각화하는 함수.
        """
        import numpy as np
        import os, glob
        import matplotlib.pyplot as plt
        from matplotlib.cm import get_cmap
        from sklearn.manifold import TSNE
    
        # [B, L, D] -> [N, D]
        ctx_np = context_codon.view(-1, context_codon.shape[-1]).detach().cpu().numpy()
        # [B, L] -> [N]
        cod_ids = seq_codon.view(-1).detach().cpu().numpy().astype(int)
    
        print(f"[codon-vis] ctx: {ctx_np.shape}, cod_ids: {cod_ids.shape}")
    
        # 1) 유효 코돈만 필터링 (22~85 / codon_id2triplet 키 기준)
        valid_ids = np.array(list(codon_id2triplet.keys()), dtype=int)
        is_valid = np.isin(cod_ids, valid_ids)
    
        if is_valid.sum() == 0:
            print("[codon-vis] No valid codon tokens to visualize.")
            return
    
        cod_ids_valid = cod_ids[is_valid]
        ctx_valid = ctx_np[is_valid]
    
        # 2) id → 코돈 문자열
        labels_np = np.array([codon_id2triplet[int(i)] for i in cod_ids_valid])
    
        print(f"[codon-vis] kept tokens: {len(labels_np)}, unique codons: {len(set(labels_np))}")
    
        # 3) 너무 많으면 샘플링
        # max_points = 8000
        # if ctx_valid.shape[0] > max_points:
        #     idx = np.random.choice(ctx_valid.shape[0], max_points, replace=False)
        #     ctx_valid = ctx_valid[idx]
        #     labels_np = labels_np[idx]





        
    
        # 4) t-SNE
        reducer = TSNE(
            n_components=2,
            perplexity=30,
            learning_rate='auto',
            init='random',
            random_state=42
        )
        reduced_np = reducer.fit_transform(ctx_valid)
    



        """
        # 4) 차원 축소: UMAP 우선, 없으면 t-SNE fallback
        try:
            import umap
            reducer = umap.UMAP(
                n_neighbors=15,     # 근접 이웃 수 (로컬 구조 강조 정도)
                min_dist=0.1,       # 점들 사이 최소 거리 (작을수록 더 촘촘히)
                n_components=2,
                metric="euclidean",
                random_state=42
            )
            print("[codon-vis] Using UMAP for dimensionality reduction.")
            reduced_np = reducer.fit_transform(ctx_valid)
        except ImportError:
            print("[codon-vis] umap-learn not installed, falling back to t-SNE.")
            from sklearn.manifold import TSNE
            reducer = TSNE(
                n_components=2,
                perplexity=30,
                learning_rate='auto',
                init='random',
                random_state=42
            )
            reduced_np = reducer.fit_transform(ctx_valid)

        """


        

        
    
        # ================= 5) 플롯 (색 = AA, 텍스트 = codon) =================
        
        
        
        # # tsne
        # os.makedirs('tsne_codon_vocab_plots_full', exist_ok=True)
        # existing = sorted(glob.glob("tsne_codon_vocab_plots_full/codon_embed_vocab_full_*.png"))
        
        
        
        
        
        # umap
        os.makedirs('umap_codon_vocab_plots_full', exist_ok=True)
        existing = sorted(glob.glob("umap_codon_vocab_plots_full/codon_embed_vocab_full_*.png"))
        
        
        
        
        
        plot_id = len(existing)
    
        fig, ax = plt.subplots(figsize=(8, 6))
    
        # (1) codon → AA 라벨 만들기
        codon_labels = labels_np                      # ex) 'GCU', 'ACC', ...
        aa_labels = np.array([codon_to_aa.get(c, '?') for c in codon_labels])
    
        # 실제로 등장한 아미노산만
        aa_set = sorted(set(aa_labels) - {'?'})
        cmap = get_cmap("tab20")
        aa2color = {aa: cmap(i % 20) for i, aa in enumerate(aa_set)}
    
        unique_codon_labels = sorted(set(codon_labels))
    
        # (2) 점 그리기: codon별로 뿌리되, 색깔은 AA 기준
        for codon in unique_codon_labels:
            mask = (codon_labels == codon)
            aa = codon_to_aa.get(codon, '?')
            if aa == '?' or aa not in aa2color:
                continue
            ax.scatter(
                reduced_np[mask, 0], reduced_np[mask, 1],
                s=8, alpha=0.4, color=aa2color[aa]
            )
    
        # (3) 텍스트 라벨: 코돈 이름을 클러스터 중앙에 찍기
        for codon in unique_codon_labels:
            mask = (codon_labels == codon)
            xs = reduced_np[mask, 0]
            ys = reduced_np[mask, 1]
            if xs.size == 0:
                continue
            cx = xs.mean()
            cy = ys.mean()
            ax.text(
                cx, cy, codon,
                fontsize=6, ha='center', va='center', weight='bold'
            )
    
        # (4) Legend: 아미노산만 표시
        handles = []
        for aa in aa_set:
            handles.append(
                plt.Line2D(
                    [0], [0],
                    marker='o',
                    color='none',
                    markerfacecolor=aa2color[aa],
                    markersize=6,
                    label=aa
                )
            )
    
        ax.set_title("Contextual Codon Embedding")
        ax.axis("off")
        ax.legend(
            handles=handles,
            title="Amino acids",
            loc='center left', bbox_to_anchor=(1.02, 0.5),
            fontsize=7, frameon=False
        )
    
        # plt.savefig(
        #     f"tsne_codon_vocab_plots_full/codon_embed_vocab_full_{plot_id:03d}.png",
        #     dpi=300, bbox_inches="tight", pad_inches=0.1
        # )
        # plt.close()
        # print(f"[✅ codon-vis] Saved codon vocab embedding plot → codon_embed_vocab_full_{plot_id:03d}.png")




        plt.savefig(
            f"umap_codon_vocab_plots_full/codon_embed_vocab_full_{plot_id:03d}.png",
            dpi=300, bbox_inches="tight", pad_inches=0.1
        )
        plt.close()
        print(f"[✅ codon-vis] Saved codon vocab embedding plot → codon_embed_vocab_full_{plot_id:03d}.png")

    



    def visualize_embedding_codon_vocab_CL(self, context_codon, seq_codon):
        """
        seq_codon (vocab ID: 22~85)을 기준으로
        codon 임베딩을 t-SNE로 시각화하는 함수.
        """
        import numpy as np
        import os, glob
        import matplotlib.pyplot as plt
        from matplotlib.cm import get_cmap
        from sklearn.manifold import TSNE
    
        # [B, L, D] -> [N, D]
        ctx_np = context_codon.view(-1, context_codon.shape[-1]).detach().cpu().numpy()
        # [B, L] -> [N]
        cod_ids = seq_codon.view(-1).detach().cpu().numpy().astype(int)
    
        print(f"[codon-vis] ctx: {ctx_np.shape}, cod_ids: {cod_ids.shape}")
    
        # 1) 유효 코돈만 필터링 (22~85 / codon_id2triplet 키 기준)
        valid_ids = np.array(list(codon_id2triplet.keys()), dtype=int)
        is_valid = np.isin(cod_ids, valid_ids)
    
        if is_valid.sum() == 0:
            print("[codon-vis] No valid codon tokens to visualize.")
            return
    
        cod_ids_valid = cod_ids[is_valid]
        ctx_valid = ctx_np[is_valid]
    
        # 2) id → 코돈 문자열
        labels_np = np.array([codon_id2triplet[int(i)] for i in cod_ids_valid])
    
        print(f"[codon-vis] kept tokens: {len(labels_np)}, unique codons: {len(set(labels_np))}")
    
        # 3) 너무 많으면 샘플링
        # max_points = 8000
        # if ctx_valid.shape[0] > max_points:
        #     idx = np.random.choice(ctx_valid.shape[0], max_points, replace=False)
        #     ctx_valid = ctx_valid[idx]
        #     labels_np = labels_np[idx]






        
        # 4) t-SNE
        reducer = TSNE(
            n_components=2,
            perplexity=30,
            learning_rate='auto',
            init='random',
            random_state=42
        )
        reduced_np = reducer.fit_transform(ctx_valid)
        




        # 4) 차원 축소: UMAP 우선, 없으면 t-SNE fallback
        """
        try:
            import umap
            reducer = umap.UMAP(
                n_neighbors=15,     # 근접 이웃 수 (로컬 구조 강조 정도)
                min_dist=0.1,       # 점들 사이 최소 거리 (작을수록 더 촘촘히)
                n_components=2,
                metric="euclidean",
                random_state=42
            )
            print("[codon-vis] Using UMAP for dimensionality reduction.")
            reduced_np = reducer.fit_transform(ctx_valid)
        except ImportError:
            print("[codon-vis] umap-learn not installed, falling back to t-SNE.")
            from sklearn.manifold import TSNE
            reducer = TSNE(
                n_components=2,
                perplexity=30,
                learning_rate='auto',
                init='random',
                random_state=42
            )
            reduced_np = reducer.fit_transform(ctx_valid)
    
        """

        
    
        # ================= 5) 플롯 (색 = AA, 텍스트 = codon) =================
        os.makedirs('tsne_codon_vocab_plots_CL_full', exist_ok=True)
        existing = sorted(glob.glob("tsne_codon_vocab_plots_CL_full/codon_embed_vocab_CL_full_*.png"))
        plot_id = len(existing)
    
        fig, ax = plt.subplots(figsize=(8, 6))
    
        # (1) codon → AA 라벨 만들기
        codon_labels = labels_np                      # ex) 'GCU', 'ACC', ...
        aa_labels = np.array([codon_to_aa.get(c, '?') for c in codon_labels])
    
        # 실제로 등장한 아미노산만
        aa_set = sorted(set(aa_labels) - {'?'})
        cmap = get_cmap("tab20")
        aa2color = {aa: cmap(i % 20) for i, aa in enumerate(aa_set)}
    
        unique_codon_labels = sorted(set(codon_labels))
    
        # (2) 점 그리기: codon별로 뿌리되, 색깔은 AA 기준
        for codon in unique_codon_labels:
            mask = (codon_labels == codon)
            aa = codon_to_aa.get(codon, '?')
            if aa == '?' or aa not in aa2color:
                continue
            ax.scatter(
                reduced_np[mask, 0], reduced_np[mask, 1],
                s=8, alpha=0.4, color=aa2color[aa]
            )
    
        # (3) 텍스트 라벨: 코돈 이름을 클러스터 중앙에 찍기
        for codon in unique_codon_labels:
            mask = (codon_labels == codon)
            xs = reduced_np[mask, 0]
            ys = reduced_np[mask, 1]
            if xs.size == 0:
                continue
            cx = xs.mean()
            cy = ys.mean()
            ax.text(
                cx, cy, codon,
                fontsize=6, ha='center', va='center', weight='bold'
            )
    
        # (4) Legend: 아미노산만 표시
        handles = []
        for aa in aa_set:
            handles.append(
                plt.Line2D(
                    [0], [0],
                    marker='o',
                    color='none',
                    markerfacecolor=aa2color[aa],
                    markersize=6,
                    label=aa
                )
            )
    
        ax.set_title("Contextual Codon Embedding (Contrastive)")
        ax.axis("off")
        ax.legend(
            handles=handles,
            title="Amino acids",
            loc='center left', bbox_to_anchor=(1.02, 0.5),
            fontsize=7, frameon=False
        )
    
        plt.savefig(
            f"tsne_codon_vocab_plots_CL_full/codon_embed_vocab_CL_full_{plot_id:03d}.png",
            dpi=300, bbox_inches="tight", pad_inches=0.1
        )
        plt.close()
        print(f"[✅ codon-vis] Saved codon vocab embedding plot → codon_embed_vocab_CL_full_{plot_id:03d}.png")


    







    def test_visualize_embedding_codon_vocab_CL(self, context_codon, seq_codon):
        """
        seq_codon (vocab ID: 22~85)을 기준으로
        codon 임베딩을 t-SNE로 시각화하는 함수.
        """
        import numpy as np
        import os, glob
        import matplotlib.pyplot as plt
        from matplotlib.cm import get_cmap
        from sklearn.manifold import TSNE
    
        # [B, L, D] -> [N, D]
        ctx_np = context_codon.view(-1, context_codon.shape[-1]).detach().cpu().numpy()
        # [B, L] -> [N]
        cod_ids = seq_codon.view(-1).detach().cpu().numpy().astype(int)
    
        print(f"[codon-vis] ctx: {ctx_np.shape}, cod_ids: {cod_ids.shape}")
    
        # 1) 유효 코돈만 필터링 (22~85 / codon_id2triplet 키 기준)
        valid_ids = np.array(list(codon_id2triplet.keys()), dtype=int)
        is_valid = np.isin(cod_ids, valid_ids)
    
        if is_valid.sum() == 0:
            print("[codon-vis] No valid codon tokens to visualize.")
            return
    
        cod_ids_valid = cod_ids[is_valid]
        ctx_valid = ctx_np[is_valid]
    
        # 2) id → 코돈 문자열
        labels_np = np.array([codon_id2triplet[int(i)] for i in cod_ids_valid])
    
        print(f"[codon-vis] kept tokens: {len(labels_np)}, unique codons: {len(set(labels_np))}")
    
        # 3) 너무 많으면 샘플링
        # max_points = 8000
        # if ctx_valid.shape[0] > max_points:
        #     idx = np.random.choice(ctx_valid.shape[0], max_points, replace=False)
        #     ctx_valid = ctx_valid[idx]
        #     labels_np = labels_np[idx]






        
        # 4) t-SNE
        reducer = TSNE(
            n_components=2,
            perplexity=30,
            learning_rate='auto',
            init='random',
            random_state=42
        )
        reduced_np = reducer.fit_transform(ctx_valid)
    
        # ================= 5) 플롯 시작 =================
        os.makedirs('test_tsne_codon_vocab_plots_CL_smalllambda', exist_ok=True)
        existing = sorted(glob.glob("test_tsne_codon_vocab_plots_CL_smalllambda/codon_embed_vocab_CL_smalllambda_*.png"))
        plot_id = len(existing)
    
        fig, ax = plt.subplots(figsize=(8, 6))
    
        codon_labels = labels_np  # 예: 'GCU', 'ACC', ...
        # codon_to_aa 는 dicts_7 같은 데서 import 되어 있다고 가정
        aa_labels = np.array([codon_to_aa.get(c, '?') for c in codon_labels])
    
        # 실제 등장한 AA만 색 부여
        aa_set = sorted(set(aa_labels) - {'?'})
        cmap = get_cmap("tab20")
        aa2color = {aa: cmap(i % 20) for i, aa in enumerate(aa_set)}
    
        unique_codon_labels = sorted(set(codon_labels))
    
        # ----- 5-1) 점 뿌리기 (색 = 아미노산) -----
        for codon in unique_codon_labels:
            mask = (codon_labels == codon)
            X = reduced_np[mask]         # [n_codon, 2]
            aa = codon_to_aa.get(codon, '?')
            if aa == '?' or aa not in aa2color:
                continue
            ax.scatter(
                X[:, 0], X[:, 1],
                s=7, alpha=0.35, color=aa2color[aa]
            )
    
        # ----- 5-2) 코돈당 여러 클러스터 center 라벨링 -----
        max_centers_per_codon = 3   # 코돈당 최대 몇 개까지 라벨 찍을지
        min_points_per_center = 80  # 이 이상 점이 있을 때만 KMeans로 여러 센터
    
        for codon in unique_codon_labels:
            mask = (codon_labels == codon)
            X = reduced_np[mask]
            if X.shape[0] == 0:
                continue
    
            n = X.shape[0]
            centers = []
    
            if n < min_points_per_center:
                # 점이 적으면 걍 평균 한 번만
                centers.append(X.mean(axis=0))
            else:
                # 코돈별 점 개수에 따라 센터 개수 결정
                k = min(max_centers_per_codon, n // min_points_per_center)
                if k <= 1:
                    centers.append(X.mean(axis=0))
                else:
                    km = KMeans(n_clusters=k, random_state=42).fit(X)
                    centers.extend(km.cluster_centers_)
    
            # 각 center에 코돈 문자열 찍기
            for cx, cy in centers:
                ax.text(
                    cx, cy, codon,
                    fontsize=6, ha='center', va='center',
                    weight='bold'
                )
    
        # ----- 5-3) legend: 아미노산 단위 -----
        handles = []
        for aa in aa_set:
            handles.append(
                plt.Line2D(
                    [0], [0],
                    marker='o',
                    color='none',
                    markerfacecolor=aa2color[aa],
                    markersize=6,
                    label=aa
                )
            )
    
        ax.set_title("Contextual Codon Embedding (Contrastive)")
        ax.axis("off")
        ax.legend(
            handles=handles,
            title="Amino acids",
            loc='center left', bbox_to_anchor=(1.02, 0.5),
            fontsize=7, frameon=False
        )
    
        plt.savefig(
            f"test_tsne_codon_vocab_plots_CL_smalllambda/codon_embed_vocab_CL_smalllambda_{plot_id:03d}.png",
            dpi=300, bbox_inches="tight", pad_inches=0.1
        )
        plt.close()
        print(f"[✅ codon-vis] Saved codon vocab embedding plot → codon_embed_vocab_CL_smalllambda_{plot_id:03d}.png")









    
    

    """
    def visualize_embedding_codon_global(self, context_codon, seq_aa, seq_codon):
        
        #t-SNE로 Codon 임베딩 시각화 (global codon dictionary 기반)
        
        
    
        # ===== 1️⃣ 기본 flatten =====
        ctx_np = context_codon.view(-1, context_codon.shape[-1]).detach().cpu().numpy()
        aa_np = seq_aa.view(-1).detach().cpu().numpy()
        cod_local_np = seq_codon.view(-1).detach().cpu().numpy()
    
        print(f"[1️⃣] context_codon: {ctx_np.shape}, seq_aa: {aa_np.shape}, seq_codon: {cod_local_np.shape}")
    
        # ===== 2️⃣ PAD 제거 =====
        PAD_AA = 21  # AA에서의 PAD 인덱스 (데이터셋에 맞게 조정)
        valid_mask = (aa_np != PAD_AA)
        aa_valid = aa_np[valid_mask]
        cod_local_valid = cod_local_np[valid_mask]
        ctx_valid = ctx_np[valid_mask]
    
        print(f"[2️⃣] after PAD filter: {len(aa_valid)} tokens left")
    
        if len(aa_valid) == 0:
            print("[⚠️] No valid tokens to visualize (maybe all PAD early in training).")
            return
    
        # ===== 3️⃣ AA 인덱스 → 문자 =====
        aa_chars = []
        miss_aa = 0
        for a in aa_valid:
            a = int(a)
            lab = index2aminoacid.get(a, None)
            if lab is None or lab == '' or lab == 'PAD':
                aa_chars.append(None)
                miss_aa += 1
            else:
                aa_chars.append(lab)
        aa_chars = np.array(aa_chars, dtype=object)
    
        print(f"[_] unique seq_aa indices (after PAD): {sorted(set(map(int, aa_valid.tolist())))[:25]} ...")
        print(f"[_] example mapping (idx->char):", [(int(a), index2aminoacid.get(int(a), None)) for a in aa_valid[:10]])
        print(f"[_] aa index not found in index2aminoacid: {miss_aa}")
    
        # ===== 4️⃣ (AA문자, local_idx) → global id =====
        global_ids = []
        miss_key = 0
        for a_char, l in zip(aa_chars, cod_local_valid):
            if a_char is None:
                global_ids.append(-1)
                continue
            key = (a_char, int(l))
            gid = aa2global_index.get(key, -1)
            if gid < 0:
                miss_key += 1
            global_ids.append(gid)
        global_ids = np.array(global_ids, dtype=int)
    
        print(f"[_] total pairs: {len(global_ids)}, unmapped(-1): {int((global_ids < 0).sum())}, miss_key count: {miss_key}")
    
        # ===== 5️⃣ 정상 매핑만 유지 =====
        keep = (global_ids >= 0)
        ctx_valid = ctx_valid[keep]
        global_ids = global_ids[keep]
    
        print(f"[_] after global mapping: {len(global_ids)} codons left")
        if len(global_ids) == 0:
            # 추가 디버그: 어떤 키들이 실패했는지 샘플 찍기
            bad = []
            for a_char, l in zip(aa_chars[:50], cod_local_valid[:50]):
                if a_char is None:
                    bad.append(("None", int(l)))
                else:
                    if (a_char, int(l)) not in aa2global_index:
                        bad.append((a_char, int(l)))
            print("[_] sample missing keys (first 20):", bad[:20])
            print("[⚠] No tokens after AA→global codon mapping. Likely index mismatch.")
            return
    
        # ===== 6️⃣ global id → codon 문자열 =====
        def to_codon_str(i: int) -> str:
            item = index2codon_global.get(int(i))
            if item is None:
                return "?"
            if isinstance(item, dict):
                return item.get("codon", "?")
            return str(item)
    
        labels_np = np.array([to_codon_str(i) for i in global_ids])
        print(f"[6️⃣] unique codon labels: {len(set(labels_np))}")
    
        # ===== 7️⃣ 샘플링 (너무 많으면 8000개만) =====
        max_points = 8000
        if ctx_valid.shape[0] > max_points:
            idx = np.random.choice(ctx_valid.shape[0], max_points, replace=False)
            ctx_valid = ctx_valid[idx]
            labels_np = labels_np[idx]
    
        # ===== 8️⃣ t-SNE =====
        reducer = TSNE(n_components=2, perplexity=30, learning_rate='auto',
                       init='random', random_state=42)
        reduced_np = reducer.fit_transform(ctx_valid)
    
        # ===== 9️⃣ 플롯 =====
        cmap = get_cmap("tab20")
        unique_labels = sorted(set(labels_np))
        label2color = {lab: cmap(i % 20) for i, lab in enumerate(unique_labels)}
    
        os.makedirs('tsne_codon_global_plots', exist_ok=True)
        existing = sorted(glob.glob("tsne_codon_global_plots/codon_embed_*.png"))
        plot_id = len(existing)
    
        fig, ax = plt.subplots(figsize=(8, 6))
        for lab in unique_labels:
            mask = (labels_np == lab)
            ax.scatter(
                reduced_np[mask, 0], reduced_np[mask, 1],
                s=10, alpha=0.6, color=label2color[lab], label=lab
            )
    
        ax.set_title("Contextual Codon Embedding (Global Mapping)")
        ax.axis("off")
        ax.legend(
            loc='center left', bbox_to_anchor=(1.02, 0.5),
            fontsize=7, markerscale=2, frameon=False
        )
    
        plt.savefig(
            f"tsne_codon_global_plots/codon_embed_{plot_id:03d}.png",
            dpi=300, bbox_inches="tight", pad_inches=0.1
        )
        plt.close()
        print(f"[✅] Saved codon global embedding plot → codon_embed_{plot_id:03d}.png")
    """







    def visualize_embedding(self, context_aa, seq_aa):
        # [B, L, D] -> [N, D]
        context_aa_np = context_aa.view(-1, context_aa.shape[-1]).detach().cpu().numpy()
        seq_aa_np = seq_aa.view(-1).detach().cpu().numpy()
    
        # 0) 디버그: raw 인덱스 분포
        unique_raw, counts_raw = np.unique(seq_aa_np, return_counts=True)
        print("[1️⃣] raw index / count")
        for idx, c in zip(unique_raw, counts_raw):
            print(f"   {idx:2d}: {c}")
    
        # 🔹 인덱스 체계 (햄이 준 버전)
        # 0~19: A..Y, 20: '*', 21: 'PAD'
        index2aminoacid = {
            0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F',
            5: 'G', 6: 'H', 7: 'I', 8: 'K', 9: 'L',
            10: 'M', 11: 'N', 12: 'P', 13: 'Q', 14: 'R',
            15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y',
            20: '*', 21: 'PAD',
        }
        PAD_ID = 21
    
        # 1) PAD만 제거 (0~20은 전부 진짜 아미노산으로 인정)
        valid_mask = (seq_aa_np != PAD_ID)
        seq_valid = seq_aa_np[valid_mask]
        ctx_valid = context_aa_np[valid_mask]
        print(f"[2️⃣] after PAD filter: {seq_valid.shape[0]} tokens left")
    
        if ctx_valid.shape[0] == 0:
            print("[⚠️] No valid amino acid tokens to visualize.")
            return
    
        # 2) 인덱스 → AA 문자로 매핑
        labels_np = np.array([index2aminoacid[i] for i in seq_valid])
    
        unique_labels, counts = np.unique(labels_np, return_counts=True)
        print("[3️⃣] label distribution after PAD filter:")
        for lab, c in zip(unique_labels, counts):
            print(f"   {lab}: {c}")
    
        # 3) 필요하면 샘플링
        max_points = 8000
        if ctx_valid.shape[0] > max_points:
            idx = np.random.choice(ctx_valid.shape[0], max_points, replace=False)
            ctx_valid = ctx_valid[idx]
            labels_np = labels_np[idx]
    
        # 4) t-SNE
        reducer = TSNE(
            n_components=2,
            perplexity=30,
            learning_rate='auto',
            init='random',
            random_state=42
        )
        reduced_np = reducer.fit_transform(ctx_valid)
        print(f"[4️⃣] reduced_np shape: {reduced_np.shape}")
    
        # 5) 플롯
        aa_order = list("ACDEFGHIKLMNPQRSTVWY*")  # 21개 (stop 포함)
        cmap = get_cmap("tab20")
        label2color = {aa: cmap(i % 20) for i, aa in enumerate(aa_order)}


        os.makedirs('highlambda_new_tsne_embedding_plots_ctx_CL_wopad', exist_ok=True)
        existing = sorted(glob.glob("highlambda_new_tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_*.png"))

        
        # os.makedirs('new_tsne_embedding_plots_ctx_CL_wopad', exist_ok=True)
        # existing = sorted(glob.glob("new_tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_*.png"))


        # # contrastive learning 없는버전 주석풀면됨
        # os.makedirs('new_tsne_embedding_plots_ctx_wopad', exist_ok=True)
        # existing = sorted(glob.glob("new_tsne_embedding_plots_ctx_wopad/context_embed_wopad_*.png"))
        
        plot_id = len(existing)
    
        fig, ax = plt.subplots(figsize=(8, 6))
    
        present_labels = sorted(set(labels_np))
        print(f"[5️⃣] Actually present AA labels: {present_labels}")
    
        for aa in aa_order:
            if aa not in present_labels:
                continue
            mask = (labels_np == aa)
            ax.scatter(
                reduced_np[mask, 0], reduced_np[mask, 1],
                color=label2color[aa],
                s=10, alpha=0.6, label=aa
            )
    
        ax.set_title('Contextual Amino Acid Embedding (CL)')
        # ax.set_title('Contextual Amino Acid Embedding')
        ax.axis('off')
        ax.legend(
            loc='center left', bbox_to_anchor=(1.02, 0.5),
            fontsize=8, markerscale=2, frameon=False
        )
    
        plt.savefig(
            f'highlambda_new_tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_{plot_id:03d}.png',
            # f'new_tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_{plot_id:03d}.png',
            # f'new_tsne_embedding_plots_ctx_wopad/context_embed_wopad_{plot_id:03d}.png',
            dpi=300, bbox_inches='tight', pad_inches=0.1
        )
        plt.close()
        print(f"[✅] Saved contextual embedding plot as highlambda_context_embed_CL_wopad_{plot_id:03d}.png")
        # print(f"[✅] Saved contextual embedding plot as context_embed_CL_wopad_{plot_id:03d}.png")
        # print(f"[✅] Saved contextual embedding plot as context_embed_wopad_{plot_id:03d}.png")







    


    
    # def visualize_embedding(self, context_aa, seq_aa):
    #     context_aa_np = context_aa.view(-1, context_aa.shape[-1]).detach().cpu().numpy()
    #     seq_aa_np = seq_aa.view(-1).detach().cpu().numpy()



        
    #     # print("aa embedding?:", context_aa_np)
    #     # print("\n embedding shape?: ", context_aa_np.shape)
    #     # print("seq_aa unique values: ", seq_aa_np)
    #     labels = [index2aminoacid[i] if i in index2aminoacid else '?' for i in seq_aa_np]
    
    #     # reducer = umap.UMAP(n_neighbors=10, min_dist=0.1, random_state=42)
    #     # reduced = reducer.fit_transform(context_aa_np)

    #     reducer = TSNE(n_components=2, perplexity=30, learning_rate='auto', init='random', random_state=42)
    #     reduced = reducer.fit_transform(context_aa_np)

    #     unique_labels = sorted(set(labels))
    #     cmap = get_cmap("tab20")
    #     label2color = {label: cmap(i % 20) for i, label in enumerate(unique_labels)}

    #     # os.makedirs('tsne_embedding_plots_ctx_wopad', exist_ok=True)
    #     os.makedirs('tsne_embedding_plots_ctx_CL_wopad', exist_ok=True)
    
    #     # ✅ 현재 폴더에 몇 개 저장돼 있는지 세기
    #     # existing = sorted(glob.glob("tsne_embedding_plots_ctx_wopad/context_embed_wopad_*.png"))
    #     existing = sorted(glob.glob("tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_*.png"))
    #     plot_id = len(existing)

        
    #     # --- plotting (replace this block) ---
    #     fig, ax = plt.subplots(figsize=(8, 6))
        
    #     labels_np = np.array(labels)
    #     reduced_np = np.array(reduced)
        
    #     # 아미노산 표준 순서(원하면 유지/수정)
    #     aa_order = list("ACDEFGHIKLMNPQRSTVWY")
    #     # 실제 존재하는 라벨만, '?'와 None 제외
    #     plot_labels = [l for l in aa_order if l in set(labels_np)]
    #     # 혹시 aa_order에 없는 라벨(예: 특수 토큰)이 있으면 뒤에 붙이기
    #     others = sorted(set(labels_np) - set(plot_labels) - set(['?']))
    #     plot_labels += others
        
    #     for label in plot_labels:
    #         if (label is None) or (label == '?'):
    #             continue
    #         mask = (labels_np == label)
    #         ax.scatter(
    #             reduced_np[mask, 0], reduced_np[mask, 1],
    #             color=label2color[label],
    #             s=10, alpha=0.6, label=label
    #         )
        
    #     ax.set_title('Contextual Amino Acid Embedding(CL)')
    #     ax.axis('off')
        
    #     # 레전드: 축 밖으로, 저장 시 잘리지 않도록 bbox 사용 예정
    #     lgd = ax.legend(
    #         loc='center left', bbox_to_anchor=(1.02, 0.5),
    #         fontsize=8, markerscale=2, frameon=False
    #     )
        
    #     # 저장 (레전드가 이미지에 포함되도록)
    #     plt.savefig(
    #         f'tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_{plot_id:03d}.png',
    #         dpi=300, bbox_inches='tight', pad_inches=0.1
    #     )
    #     plt.close()
    #     print(f"[✅] Saved contextual embedding plot as tsne_context_embed_CL_wopad_{plot_id:03d}.png")

    
        # 3. 아미노산 고정 색상 매핑
        # unique_labels = sorted(set(labels))
        # cmap = get_cmap("tab20")
        # label2color = {label: cmap(i % 20) for i, label in enumerate(unique_labels)}
        
        # # os.makedirs('tsne_embedding_plots_ctx_wopad', exist_ok=True)
        # os.makedirs('tsne_embedding_plots_ctx_CL_wopad', exist_ok=True)
    
        # # ✅ 현재 폴더에 몇 개 저장돼 있는지 세기
        # # existing = sorted(glob.glob("tsne_embedding_plots_ctx_wopad/context_embed_wopad_*.png"))
        # existing = sorted(glob.glob("tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_*.png"))
        # plot_id = len(existing)
        
        # plt.figure(figsize=(8, 6))
        # for i, label in enumerate(labels):
        #     if label == '?' or label is None:  # 라벨 없는 경우
        #         continue
        #     color = label2color[label]
        #     mask = [l == label for l in labels]
        #     plt.scatter(
        #         np.array(reduced)[mask, 0],
        #         np.array(reduced)[mask, 1],
        #         color=color,
        #         s=10,
        #         alpha=0.6,
        #         label=label  # ✅ legend용 라벨
        #     )
        #     # plt.scatter(reduced[i, 0], reduced[i, 1], color=color, s=10, alpha=0.6)
            
        # # plt.title('Contextual Amino Acid Embedding')
        # plt.title('Contextual Amino Acid Embedding(CL)')
        # plt.axis('off')
        # plt.tight_layout()

        # # legend 추가
        # plt.legend(
        #     loc='center left',
        #     bbox_to_anchor=(1.05, 0.5),  # plot 밖에 배치
        #     fontsize=8,
        #     markerscale=2  # legend 점 크기 키움
        # )
            
        # # # ✅ 자동 넘버링된 파일명으로 저장
        # # plt.savefig(f'tsne_embedding_plots_ctx_wopad/context_embed_wopad_{plot_id:03d}.png', dpi=300)
        # # plt.close()
        # # print(f"[✅] Saved contextual embedding plot as tsne_context_embed_wopad_{plot_id:03d}.png")

        # # ✅ 자동 넘버링된 파일명으로 저장
        # plt.savefig(f'tsne_embedding_plots_ctx_CL_wopad/context_embed_CL_wopad_{plot_id:03d}.png', dpi=300)
        # plt.close()
        # print(f"[✅] Saved contextual embedding plot as tsne_context_embed_CL_wopad_{plot_id:03d}.png")

    
    def step_epoch(self, data_loader, islast=False):
        losses_epoch = [0.0, 0.0, 0.0, 0.0, 0.0]
        score_epoch = [0.0, 0.0, 0.0, 0.0, 0.0]

        for loader in data_loader:
            losses, _, _, labels, preds, pads = self.step_loader(loader, islast)
            score = self.get_scores(labels, preds, pads)
            
            losses_epoch = [losses_epoch[i] + losses[i] for i in range(len(losses))]
            score_epoch = [score_epoch[i] + score[i] for i in range(len(score))]

        losses_epoch = [losses_epoch[i]/len(data_loader) for i in range(len(losses_epoch))]
        score_epoch = [score_epoch[i]/len(data_loader) for i in range(len(score_epoch))]
        return losses, score_epoch
        
    # def step_loader(self, loader, islast=False):
    #     seqs_aa, seqs_codon, seqs_label, pads, seqs_pred = torch.empty(0), torch.empty(0), torch.empty(0), torch.empty(0), torch.empty(0)
    #     losses = [0.0, 0.0, 0.0, 0.0, 0.0]  # Total, AA, Codon, Contrastive, mfe
    #     for i, data in enumerate(loader):
    #         loss, loss_aa, loss_codon, loss_contrastive, loss_mfe = 0.0, 0.0, 0.0, 0.0, 0.0


    #         seq_aa, seq_codon, seq_label, pad = data[0].to(self.device), data[1].to(self.device), data[2].to(self.device), data[3].to(self.device)
    #         if self.train_mode:
    #             self.optim.zero_grad()
            
    #         seq_pred, context_aa = self.forward(seq_aa)
    #         loss_aa = self.get_loss(seq_pred, seq_label, pad)
    #         loss += loss_aa


    #         # ✅ 항상 codon embedding 계산 (codon embedding 뽑기 위해)
    #         emb_codon = self.model.embedding(seq_codon)
    #         emb_codon = self.model.pe(emb_codon)
    #         emb_out, context_codon = self.model.backbone(emb_codon)
                        

    #         if self.train_mode:
    #             lambda_value = self.lambda_codon + self.lambda_emb + self.lambda_mfe
    #             if lambda_value > 0.0:

                    
    #                 # loss using codon sequence
    #                 if self.lambda_codon + self.lambda_emb > 0.0:
    #                     emb_codon = self.model.embedding(seq_codon)
    #                     emb_codon = self.model.pe(emb_codon)
    #                     emb_out, context_codon = self.model.backbone(emb_codon)

    #                     if self.lambda_codon > 0.0:
    #                         seq_pred_codon = self.model.to_out(emb_out)
    #                         loss_codon = self.get_loss(seq_pred_codon, seq_label, pad)
    #                         loss = loss + (self.lambda_codon * loss_codon)

    #                     if self.lambda_emb > 0.0:
    #                         context_aa = self.model.to_proj(context_aa)
    #                         context_codon = self.model.to_proj(context_codon)
    #                         loss_contrastive = self.CELOSS_window(context_aa, context_codon, tau=0.5, window=self.window, step_size=self.step_size)
    #                         loss = loss + (self.lambda_emb * loss_contrastive)

    #                 # loss using MFE
    #                 if self.lambda_mfe > 0.0:
    #                     loss_mfe = self.get_loss_mfe(seq_pred, seq_label, seq_aa, pad, islast)
    #                     loss = loss + (self.lambda_mfe * loss_mfe)
                
    #             loss.backward()
    #             self.optim.step()
            
    #         losses[0] += loss
    #         losses[1] += loss_aa
    #         losses[2] += loss_codon
    #         losses[3] += loss_contrastive
    #         losses[4] += loss_mfe
            
    #         if self.verbose:
    #             print('\r', self.ticks[(i%4)], f'[{i+1}/{len(loader)}] Loss: {loss:.5f} [batch: {len(seq_aa)} length: {len(seq_aa[0])}]', end = '')

    #         seqs_aa = torch.cat([seqs_aa, seq_aa.cpu()], dim=0)
    #         seqs_codon = torch.cat([seqs_codon, seq_aa.cpu()], dim=0)
    #         seqs_label = torch.cat([seqs_label, seq_label.cpu()], dim=0)
    #         seqs_pred = torch.cat([seqs_pred, seq_pred.cpu()], dim=0)
    #         pads = torch.cat([pads, pad.cpu()], dim=0)          







        
    #     # 여기서 embedding visualize하기 (for문 밖에서 해야함)
    #     if islast:
    #         # self.visualize_embedding(context_aa, seq_aa)

    #         self.visualize_embedding_codon(
    #         context_codon=context_codon,        # [B, L, D]
    #         seq_codon=seq_codon,                # [B, L]
    #         index2codon=self.data.index2codon,  # 또는 tokenizer.index2codon
    #         pad_id_codon=self.data.pad_id_codon,
    #         title="Contextual Codon Embedding (CL)" if self.lambda_emb>0 else "Contextual Codon Embedding",
    #         out_dir="tsne_codon"  # 폴더명 취향대로
    #     )
        
        
        
    #     losses = [losses[i]/len(loader) for i in range(len(losses))]
    #     return losses, seqs_aa, seqs_codon, seqs_label, seqs_pred, pads




















    def step_loader(self, loader, islast=False):
        seqs_aa, seqs_codon, seqs_label, pads, seqs_pred = torch.empty(0), torch.empty(0), torch.empty(0), torch.empty(0), torch.empty(0)
        losses = [0.0, 0.0, 0.0, 0.0, 0.0]
    
        # 시각화용 누적 버퍼 (마지막에 한 번에 그림)
        all_ctx_codon, all_seq_aa, all_seq_codon = [], [], []
    
        for i, data in enumerate(loader):
            loss = loss_aa = loss_codon = loss_contrastive = loss_mfe = 0.0
    
            seq_aa   = data[0].to(self.device)     # [B, L]  (AA 인덱스)
            seq_codon= data[1].to(self.device)     # [B, L]  (AA-내 local codon 인덱스)
            seq_label= data[2].to(self.device)
            pad      = data[3].to(self.device)
    
            if self.train_mode:
                self.optim.zero_grad()
    
            # AA 스트림
            seq_pred, context_aa = self.forward(seq_aa)
            loss_aa = self.get_loss(seq_pred, seq_label, pad)
            loss += loss_aa
    
            # ✅ codon 컨텍스트는 항상 계산 (loss에만 조건)
            emb_codon = self.model.embedding(seq_codon)
            emb_codon = self.model.pe(emb_codon)
            emb_out, context_codon = self.model.backbone(emb_codon)
    
            if self.train_mode:
                # codon loss
                if self.lambda_codon > 0.0:
                    seq_pred_codon = self.model.to_out(emb_out)
                    loss_codon = self.get_loss(seq_pred_codon, seq_label, pad)
                    loss += self.lambda_codon * loss_codon
    
                # contrastive loss
                if self.lambda_emb > 0.0:
                    ctx_aa_proj    = self.model.to_proj(context_aa)
                    ctx_codon_proj = self.model.to_proj(context_codon)
                    loss_contrastive = self.CELOSS_window(
                        ctx_aa_proj, ctx_codon_proj, tau=0.5,
                        window=self.window, step_size=self.step_size
                    )
                    loss += self.lambda_emb * loss_contrastive
    
                # MFE
                if self.lambda_mfe > 0.0:
                    loss_mfe = self.get_loss_mfe(seq_pred, seq_label, seq_aa, pad, islast)
                    loss += self.lambda_mfe * loss_mfe
    
                loss.backward()
                self.optim.step()
    
            # 로그/누적
            losses[0] += loss; losses[1] += loss_aa; losses[2] += loss_codon
            losses[3] += loss_contrastive; losses[4] += loss_mfe
    
            if self.verbose:
                print('\r', self.ticks[(i%4)], f'[{i+1}/{len(loader)}] Loss: {loss:.5f} [batch: {len(seq_aa)} length: {len(seq_aa[0])}]', end='')
    
            seqs_aa   = torch.cat([seqs_aa,   seq_aa.cpu()], dim=0)
            seqs_codon= torch.cat([seqs_codon,seq_codon.cpu()], dim=0)
            seqs_label= torch.cat([seqs_label,seq_label.cpu()], dim=0)
            seqs_pred = torch.cat([seqs_pred, seq_pred.cpu()], dim=0)
            pads      = torch.cat([pads,      pad.cpu()], dim=0)
    
            # ✅ 시각화용 누적
            all_ctx_codon.append(context_codon.detach().cpu())
            all_seq_aa.append(seq_aa.detach().cpu())
            all_seq_codon.append(seq_codon.detach().cpu())
    
        # ✅ for문 끝난 뒤 딱 한 번 시각화
        # if islast and context_codon is not None:
        #     # self.visualize_embedding_codon_vocab(context_codon, seq_codon)
        #     # self.visualize_embedding_codon_vocab_CL(context_codon, seq_codon)
        #     self.test_visualize_embedding_codon_vocab_CL(context_codon, seq_codon)




        











        
        return losses, seqs_aa, seqs_codon, seqs_label, seqs_pred, pads


















    
    def get_loss_mfe(self, preds, labels, seq_amino, n_pads, islast=False):
        seqs_pred, seqs_amino = [], []
        values_mfe, values_cai = [], []

        # Apply argmax before detaching to ensure integer sequence prediction
        preds_argmax = torch.argmax(preds, dim=-1)

        # Move tensors to CPU and detach for better efficiency
        preds_np = preds_argmax.cpu().detach().numpy()
        seq_amino_np = seq_amino.cpu().detach().numpy()
        n_pads_np = n_pads.cpu().numpy()

        if self.wo_pad:
            # Adjust predictions and amino acid sequences by removing padding
            for i, n_pad in enumerate(n_pads_np):
                if n_pad > 0:
                    seq_pred_trimmed = preds_np[i, :-int(n_pad)]
                    seq_amino_trimmed = seq_amino_np[i, :-int(n_pad)]
                else:
                    seq_pred_trimmed = preds_np[i, :]
                    seq_amino_trimmed = seq_amino_np[i, :]

                seqs_pred.append(seq_pred_trimmed)
                seqs_amino.append(seq_amino_trimmed)

                # Calculate MFE and CAI per sequence
                mfe_tmp, cai_tmp = self.optimizer._calculate_loss_mfe_cai(seq_amino_trimmed, seq_pred_trimmed)
                values_mfe.append(mfe_tmp)
                values_cai.append(cai_tmp)
        else:
            seqs_pred = list(preds_np)
            seqs_amino = list(seq_amino_np)

            # Vectorized calculation for all sequences at once (if supported)
            values_mfe, values_cai = zip(*[
                self.optimizer._calculate_loss_mfe_cai(seq_amino_np[i], preds_np[i]) 
                for i in range(len(seq_amino_np))
            ])

        # Handle empty lists to prevent division by zero
        if len(values_mfe) > 0:
            mfe = np.mean(values_mfe)
            cai = np.mean(values_cai)
        else:
            mfe, cai = 0.0, 0.0

        # Convert to PyTorch tensors
        mfe_tensor = torch.tensor(mfe, dtype=torch.float32, device=self.device)
        cai_tensor = torch.tensor(cai, dtype=torch.float32, device=self.device)

        # # Debugging information
        # if self.verbose:
        #     print(f"[DEBUG] MFE Loss (per seq): {values_mfe}")
        #     print(f"[DEBUG] CAI Loss (per seq): {values_cai}")
        #     print(f"[DEBUG] MFE Loss (avg): {mfe}, CAI Loss (avg): {cai}")

        # if islast and self.verbose:
        #     print(f"[DEBUG] Final MFE Loss: {mfe_tensor.item()}, Final CAI Loss: {cai_tensor.item()}")
        weight_mfe = 4.0
        return torch.sigmoid((weight_mfe* mfe_tensor - cai_tensor)/2)

    def forward(self, x):
        n_iter=1
        if self.iter_mcdropout > 0:
            n_iter = 1 if self.train_mode else self.iter_mcdropout

        emb = self.model.embedding(x)
        emb = self.model.pe(emb)
        emb_out, emb_context = self.model.backbone(emb)
        out = self.model.to_out(self.mcdropout(emb_out))

        for _ in range(n_iter-1):
            out += self.model.to_out(self.mcdropout(emb_out))

        return out / n_iter, emb_context
    
    def without_pads(self, pred, codon, n_pads, amino=None):
        amino_wo_pads, pred_wo_pads, codon_wo_pads = [], [], []
        for i, n_pad in enumerate(n_pads):
            if n_pad > 0:
                if amino is not None:
                    amino_wo_pads.append(amino[i, :-int(n_pad.item())])
                codon_wo_pads.append(codon[i, :-int(n_pad.item())])
                pred_wo_pads.append(pred[i, :-int(n_pad.item())])
            else:
                if amino is not None:
                    amino_wo_pads.append(amino[i, :])
                pred_wo_pads.append(pred[i, :])
                codon_wo_pads.append(codon[i, :])

        if amino is not None:
            amino_wo_pads = torch.cat(amino_wo_pads).to(self.device)
        else:
            amino_wo_pads = None

        pred_wo_pads = torch.cat(pred_wo_pads).to(self.device)
        codon_wo_pads = torch.cat(codon_wo_pads).to(self.device)
        return pred_wo_pads, codon_wo_pads, amino_wo_pads

    def get_loss(self, predictions, targets, n_pads):
        if self.wo_pad:
            predictions, targets, _ = self.without_pads(predictions, targets, n_pads)
        else:
            _, _, d = predictions.shape
            predictions = predictions.view(-1, d)
            targets = targets.view(-1)

        loss = self.loss_fn(predictions, targets)
        return loss

    def get_scores(self, codons, preds, n_pads):
        print('\r Calculating...', end='')
        if self.wo_pad:
            preds, codons, _ = self.without_pads(preds, codons, n_pads)
        
        print('ACC ', end=' ')
        codons_np = codons.view(-1).long().cpu().detach().numpy()
        preds_np = torch.argmax(preds, dim=-1).view(-1).cpu().numpy()
        correct = preds_np == codons_np
        acc = correct.mean()

        print('f1m ', end=' ')
        f1_ma = f1_score(codons_np, preds_np, average='macro')
        print('f1w ', end=' ')
        f1_w = f1_score(codons_np, preds_np, average='weighted')
        print('MCC ', end=' ')
        mcc = matthews_corrcoef(codons_np, preds_np)
        return [acc, f1_ma, f1_w, mcc]

    def get_statistics(self, df):
        df = df[['Name', 'AminoAcid']].rename(columns={'Name':'Name', 'AminoAcid':'Sequence'})
        df['Type'] = ['Protein']*df.shape[0]
        df = df[['Name', 'Sequence', 'Type']]
        df = df.drop_duplicates(subset='Name')
        df = df.reset_index(drop=True)

        print('\r Calculating...', end='')
        df_stat = self.optimizer.get_optimized_codons(self.model, df)
        df_score = self.optimizer._calculate_codon_optimization_metrics(df_stat)
        print('CAI ',end=' ')
        cai=df_score['CAI'].mean()
        print('GC ',end=' ')
        gc=df_score['GCContent'].mean()
        print('GC3 ',end=' ')
        gc3=df_score['GC3'].mean()
        #print('MFE ',end=' ')
        #mfe=df_score['MFE'].mean()

        #return [cai, gc, gc3, mfe]
        return [cai, gc, gc3]
    
    def get_hgcOPT_statistics(self, df):
        df = df[['Name', 'AminoAcid']].rename(columns={'Name':'Name', 'AminoAcid':'Sequence'})
        df['Type'] = ['Protein']*df.shape[0]
        df = df[['Name', 'Sequence', 'Type']]
        df = df.drop_duplicates(subset='Name')
        df = df.reset_index(drop=True)

        print('\r Calculating...', end='')
        df_stat = self.optimizer.get_optimized_codons(self.model, df)
        df_score = self.optimizer._calculate_hgcOPT(df_stat)
        rna = df_score['OptimizedRNA'].values[0]
        codon = df_score['OptimizedCodon'].values[0]
        print('CAI ',end=' ')
        cai=df_score['CAI'].values[0]
        print('GC ',end=' ')
        gc=df_score['GCContent'].values[0]
        print('GC3 ',end=' ')
        gc3=df_score['GC3'].values[0]
        
        print('MFE ',end=' ')
        mfe=df_score['MFE'].values[0]

        return [cai, gc, gc3, mfe, rna, codon]
        # return [cai, gc, gc3, rna, codon]

    def CELOSS_window(self, out1, out2, tau=0.5, window=16, step_size=16):
        total_loss = 0.0
        num_windows = 0

        seq_len = out1.size(1)
        for i in range(0, seq_len, step_size):
            if i + window <= seq_len:
                out1_window = out1[:, i:i + window, :]
                out2_window = out2[:, i:i + window, :]
            else:
                # Handle the last window if seq_len is not perfectly divisible by window size
                out1_window = out1[:, -window:, :]
                out2_window = out2[:, -window:, :]

            loss = self.CELoss(out1_window, out2_window, tau)
            total_loss += loss
            num_windows += 1

            # Debugging output
            if torch.any(loss < 0):
                print(f"Negative loss detected in window [{i}:{i + window}]. Loss: {loss.item()}")

        # Ensure the total loss is non-negative
        final_loss = total_loss / num_windows if num_windows > 0 else 0.0
        if final_loss < 0:
            print("Final loss is negative. This should not happen.")
            print("total_loss:", total_loss)
            print("num_windows:", num_windows)
            print("final_loss:", final_loss)

        return final_loss
    
    def CELoss(self, out1, out2, tau=0.5):
        eps = 1e-10
        sim = F.cosine_similarity(out1, out2, dim=2)
        
        # Add numerical stability check for sim values
        sim = torch.clamp(sim, min=-1 + eps, max=1 - eps)

        exp_sim = torch.exp(sim / tau)

        # Ensure numerical stability
        exp_sim_sum = exp_sim.sum(dim=1, keepdim=True) + eps
        exp_sim_diagonal = torch.diagonal(exp_sim, dim1=-2, dim2=-1) + eps

        # Compute contrasts and handle potential negative values in logarithm
        contrasts = -torch.log(exp_sim_diagonal / exp_sim_sum)
        contrasts = torch.clamp(contrasts, min=0.0)

        # Compute the mean loss
        loss = contrasts.mean()
        return loss