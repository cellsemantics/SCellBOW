from opcode import stack_effect
from imblearn.over_sampling import SMOTE
from SCellBOW_cluster import SCellBOW_clust
from sksurv.ensemble import RandomSurvivalForest
from collections import defaultdict
from sklearn.preprocessing import MinMaxScaler
from nltk.tokenize import word_tokenize
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import time
import os
import random
from tqdm import tqdm
import matplotlib.pyplot as plt
import nltk
from xgbse.converters import convert_to_structured
nltk.download('punkt')
os.environ['PYTHONHASHSEED'] = '0'  # addition


# class SCellBOW_algebra(SCellBOW_clust):
#     def __init__(self, data_target, save_dir, iter=40):

#         SCellBOW_clust.__init__(self, adata_target, save_dir, iter=40)

#     def get_from_Cluster(self):
#         p = self.SCellBOW_target()
#         return p




class SCellBOW_algebra(SCellBOW_clust):
    # Transfer Learning
    def __init__(self, adata_test, adata_train, save_dir, Type='clusters', iter=40, bootstrap_samples=50, split=0.2, unit="UMI", n_top_features=1000):
        self.adata_test = adata_test  # Load single cell data
        self.adata_train = adata_train  # Load Bulk RNAseq data
        self.save_dir = "./"+save_dir+"/"
        if not os.path.exists(self.save_dir):
            raise Exception("Train Model not Found")
        self.iter = iter
        self.big_model = Doc2Vec.load(self.save_dir + 'source_model')
        if not os.path.isfile(self.save_dir+'source_model'):
            raise Exception("Source Model not found")
        self.Type = Type
        self.bootstrap_samples = bootstrap_samples
        self.split = split
        self.unit = unit
        self.n_top_features = n_top_features
        self.SCellBOW_phenotype_algebra()
        
    def load(self, path):
        dbfile = open(path, 'rb')
        data = pickle.load(dbfile)
        dbfile.close()
        return data

    def save(self,data, path):
        dbfile = open(path, 'wb')
        pickle.dump(data, dbfile, protocol=4)
        dbfile.close()
        return None

    def get_from_SCellBOW_cluster(self):
        p = self.SCellBOW_target()
        return p
    
    

    # Create document out of Gene Expression matrix

    # Shuffle corpus

    def class_imbalance(self, adata, col):
        matrix = adata.to_df()
        matrix["type"] = adata.obs[col].values
        oversample = SMOTE(random_state=42)
        X, y = oversample.fit_resample(
            matrix.iloc[:, :-2], matrix['type'].values)
        unique, counts = np.unique(y, return_counts=True)
        dataframe = pd.DataFrame(X)
        dataframe['type'] = y
        return dataframe

    def sc_pseudobulk_anndata(self, data):
        # pseudobulk for combined type
        all_vec = data.mean(axis=0, numeric_only=True)
        all_vec = pd.DataFrame(all_vec).T
        all_vec.rename(index={0: 'All'}, inplace=True)

        # pseudobulk by individual type
        celltype_avg_vec = data.groupby('type').mean()
        cell_type = celltype_avg_vec.index  # save the cell types
        print("Shape of pseudobulk by individual type:", celltype_avg_vec.shape)

        # Add combined

        # stack All and type wise Average Vector
        stacked_vec = pd.concat([all_vec, celltype_avg_vec])

        # Create scanpy object
        adata_stack = sc.AnnData(stacked_vec)
        adata_stack.obs['type'] = stacked_vec.index.values
        return adata_stack

    def preprocessing(self, adata, unit, top_gene_num):
        adata.var_names_make_unique()  # takes unique genes
        sc.pp.filter_genes(adata, min_cells=20)
        if unit == "UMI" or unit == "CPM":
            sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        print("pre:", adata.shape)
        sc.pp.highly_variable_genes(adata, n_top_genes=top_gene_num)
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]
        print("post:", adata.shape)
        sc.pp.scale(adata, max_value=10)
        return adata

    def SCellBOW_phenotype_algebra(self):
        
      # Class balance of Single cell based on type
        data = self.class_imbalance(self.adata_test, self.Type)

        # Prepare Pseudobulk from single cell
        adata_SCell = self.sc_pseudobulk_anndata(data)

        # Combine Bulk and single-cell dataset
        adata_target = self.adata_train.concatenate(
            adata_SCell, join='inner', fill_value=0)

        # Preprocess combined dataset
        adata_target = self.preprocessing(
            adata_target, self.unit, self.n_top_features)

        # Transfer Learning
        # Call ScellBOW_test
        # INherit ScellBOW_clust
        SCellBOW_clust.__init__(self, adata_target, self.save_dir,
                                      iter=self.iter)

        adata_tl = self.get_from_SCellBOW_cluster()
        
        #self.SCellBOW_target(adata_target, self.save_dir,
        #                              iter=self.iter)
        
        print(adata_tl)
        # Load the embedding
        vector = adata_tl.obsm['SCellBOW_embed']

        # length of bulk data
        n = len(self.adata_train.to_df())

        # bulk vector
        bulk_vec = vector[:n]
        bulk_obj = sc.AnnData(bulk_vec)
        bulk_obj.obs = self.adata_train.obs
        print(bulk_obj.shape)

        # single cell vector
        single_vec = vector[n:]
        single_obj = sc.AnnData(single_vec)
        single_obj.obs = adata_SCell.obs
        print(single_obj.shape)

        # Survival Analysis
        bulk_data = pd.DataFrame(bulk_vec)
        bulk_data['duration'] = self.adata_train.obs.time.values
        bulk_data['event'] = self.adata_train.obs.status.values
        print(bulk_data.head(5))

        # Separate Single-cell data for Testing

        single_data = single_obj.to_df()
        single_data['type'] = adata_SCell.obs['type'].values
        print("Shape:", single_data.shape)
        print(single_data.head(5))

        # Prepare All celltype vs Rest celltype data
        all_vec = single_data.loc[single_data['type'] == 'All']
        all_vec = all_vec.drop(['type'], axis=1)
        rest_vec = single_data.drop(all_vec.index.values)

        # Survival Analysis
        # to store the predicted score for each patient
        predicted_score = defaultdict(list)
        # Bootstrap on rest of the data
        for j in range(0, self.bootstrap_samples):  # run boostrapping
            df_train = bulk_data
            print("loop ", j)

            df_test = df_train.sample(frac=self.split, random_state=j)
            df_train = df_train.drop(df_test.index)

            X_train = df_train.drop(['duration', 'event'], axis=1)
            y_train = convert_to_structured(
                df_train['duration'], df_train['event'])
            print(X_train.shape, y_train.shape)

            # Train the model
            pr = RandomSurvivalForest(
                random_state=0, n_jobs=-1).fit(X_train, y_train)

            # Score for combined pseudobulk
            X_test = pd.DataFrame(all_vec.values)
            key = "Combined"
            y_predicted = pr.predict(X_test)
            predicted_score[key].append(y_predicted.item(0))

            # Prediction
            for index, vector in rest_vec.iterrows():
                key = "Combined - (" + vector['type']+")"
                test_celltype = pd.DataFrame(vector.drop(['type']).values).T
                sub_test = pd.DataFrame(all_vec.values - test_celltype.values)
                predicted_score[key].append(pr.predict(sub_test).item(0))

        predicted_risk_score = pd.DataFrame(predicted_score)
        # save the predicted risk score
        return predicted_risk_score



