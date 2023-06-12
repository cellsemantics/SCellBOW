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


   

class SCellBOW_pretrain:

    def __init__(self, adata_source, save_dir, vec_size=300, n_worker=1, iter=20):
        self.adata_source = adata_source
        self.corpus_not_shuffled_trn = None
        self.corpus_shuffled_trn = None

        # self.save_dir = os.makedirs("./"+save_dir+"/") if not "./"+save_dir+"/" else "./"+save_dir+"/"
        self.save_dir = os.makedirs(
            './'+save_dir+'/') if not os.path.exists('./'+save_dir+'/') else './'+save_dir+'/'
        self.vec_size = vec_size
        self.n_worker = n_worker
        self.iter = iter
        self.run(self.iter, self.vec_size)
        

    def save(self, data, path):
        dbfile = open(path, 'wb')
        pickle.dump(data, dbfile, protocol=4)
        dbfile.close()
        return None


# Create corpus

    def wordfreq_docs_list(self, df):
        corpus = []
        row = []
        s = ''
        names = {e: name for e, name in enumerate(df.columns)}
        for i in tqdm(iterable=df.itertuples()):
            for e, j in enumerate(i[1:]):
                temp = names[e]
                row += [temp] * int(j)
            corpus.append(row)
            s = ''
            row = []
        print('corpus created with size: ', len(corpus))
        return corpus


#Shuffle corpus

    def shuf(self, cnsl):
        corpus_shuffled = []
        random.seed(0)
        for l in tqdm(range(len(cnsl))):
            random.shuffle(cnsl[l])
            s = ' '.join(cnsl[l])
            corpus_shuffled.append(s)
        return corpus_shuffled

    #Train Doc2vec


    def doc2vec(self, corpus, iter, model_name, vec_size):
        # tagging docs
        print('tagging docsX')
        tagged_data = [TaggedDocument(words=word_tokenize(
            _d), tags=[str(i)]) for i, _d in enumerate(corpus)]
        print('all docs tagged with len', len(tagged_data))

        model = Doc2Vec(vector_size=vec_size,
                  alpha=0.025,  # initial learning rate
                  min_alpha=0.00025,
                  min_count=1,
                  window=2,
                  workers=self.n_worker,
                  seed=0,
                  dm=1)

        # Building a vocabulary
        model.build_vocab(tagged_data, update=False)
        print('vocab built')
        model.train(tagged_data,
                    total_examples=model.corpus_count,
                    epochs=iter)
        # Save Model
        model.save(self.save_dir + model_name)
        print('\nmodel trained')  # add loggings
        return None
   

    def run(self, iter, vec_size):
        
        print("The save directory is" , self.save_dir)
        
        #rescale the data
        srcdata = self.adata_source.to_df()
        print(srcdata.shape)  # add to logging add timestamp

        scaler = MinMaxScaler(feature_range=(1, 10))
        print(scaler.fit(srcdata))
        trainData = scaler.transform(srcdata)
        trainData = pd.DataFrame(trainData, columns=srcdata.columns)
        # Shuffle the corpus
        self.corpus_not_shuffled_trn = self.wordfreq_docs_list(trainData)
        print('Startb Corpus shuffle')
        self.corpus_shuffled_trn = self.shuf(self.corpus_not_shuffled_trn)
        print('Start Training Model')
        #name = self.save_dir + 'source_corpus.pickle'
        # self.save(self.corpus_shuffled_trn,name )  Need work
        # Train the model
        self.doc2vec(self.corpus_shuffled_trn, iter=self.iter,
                    model_name='source_model', vec_size=vec_size)
        return None
    
    

    

class SCellBOW_cluster:
    #Transfer Learning
    def __init__(self,adata_target,save_dir,iter=40,resolution=1.0,neighbors=15):
        self.adata_target = adata_target
        self.adata = None
        self.save_dir = "./"+save_dir+"/"
        if not os.path.exists(self.save_dir):
            raise Exception("Train Model not Found")
        self.iter = iter
        self.big_model = Doc2Vec.load(self.save_dir + 'source_model')
        if not os.path.isfile(self.save_dir+'source_model'):
            raise Exception("Source Model not found")
        self.resolution = resolution
        self.neighbors = neighbors
        
    
    def save(self, data, path):
        dbfile = open(path, 'wb')
        pickle.dump(data, dbfile, protocol=4)
        dbfile.close()
        print(path)
        return None

    def wordfreq_docs_list(self, df):
        corpus = []
        row = []
        s = ''
        names = {e: name for e, name in enumerate(df.columns)}
        for i in tqdm(iterable=df.itertuples()):
            for e, j in enumerate(i[1:]):
                temp = names[e]
                row += [temp] * int(j)
            corpus.append(row)
            s = ''
            row = []
        print('corpus created with size: ', len(corpus))
        return corpus

    def shuf(self, cnsl):
        corpus_shuffled = []
        random.seed(0)
        for l in tqdm(range(len(cnsl))):
            random.shuffle(cnsl[l])
            s = ' '.join(cnsl[l])
            corpus_shuffled.append(s)
        return corpus_shuffled

    def run(self):
        print("SCellBOW_target")
        self.adata = self.adata_target.copy()
        dstdata = self.adata.to_df()
        print(dstdata.shape)
        scaler = MinMaxScaler(feature_range=(1, 10))
        print(scaler.fit(dstdata))
        trainData = scaler.transform(dstdata)
        trainData = pd.DataFrame(trainData, columns=dstdata.columns)
        print(trainData.shape)

        corpus_not_shuffled = self.wordfreq_docs_list(trainData)
        print('Corpus Not Shuffled', len(corpus_not_shuffled))
        corpus_shuffled = self.shuf(corpus_not_shuffled)
        print('Corpus Shuffled', len(corpus_shuffled))
        # tokenizing the corpus
        tagged_data_tl = [TaggedDocument(words=word_tokenize(
            _d), tags=[str(i)]) for i, _d in enumerate(corpus_shuffled)]
        print('all docs tagged with len', len(tagged_data_tl))
        self.save(tagged_data_tl, self.save_dir + 'tagged_data')

        # Updating the vocabulary
        self.big_model.build_vocab(tagged_data_tl, progress_per=1000, update=True)
        print('vocab updated')

        # Retraining the model

        self.big_model.train(
            tagged_data_tl, total_examples=self.big_model.corpus_count, epochs=self.iter)

        # infer vectors for new corpus
        re_vecs_W_new = []
        for c in tqdm(corpus_shuffled):
            c = c.split(' ')
            self.big_model.random.seed(0)  # addition
            re_vecs_W_new.append(self.big_model.infer_vector(
                c, epochs=self.iter))  # addition

        re_vecs_W_new = np.array(re_vecs_W_new)

        self.save(re_vecs_W_new, self.save_dir+'vector_embeddings')
        self.adata.obsm['SCellBOW_embed'] = re_vecs_W_new
        self.adata.obsm['X_embed'] = sc.tl.pca(
            re_vecs_W_new, svd_solver='arpack')

        sc.pp.neighbors(self.adata, n_neighbors=self.neighbors, use_rep = 'X_embed', random_state=0)
        print(self.adata)
        sc.tl.umap(self.adata)

        sc.tl.leiden(self.adata, key_added='clusters_' +str(self.resolution), resolution=self.resolution)
        print("done leiden")

        self.adata.write(self.save_dir+"adata_scellbow.h5ad")
        # with plt.rc_context({'figure.figsize': (5, 5)}):
            # sc.pl.umap(self.adata_target, color='clusters_'+str(self.resolution),show=True)
        print(self.adata)
        #adata = self.adata_target.copy()
        return self.adata



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

    def run(self):
        
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
        median_risk_score = predicted_risk_score.median()
        # save the predicted risk score
        return median_risk_score, predicted_risk_score
    
    
 
