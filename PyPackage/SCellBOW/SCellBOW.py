from opcode import stack_effect
from imblearn.over_sampling import SMOTE
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
#from tqdm import tqdm
from tqdm.notebook import tqdm_notebook
import matplotlib.pyplot as plt
import nltk
from xgbse.converters import convert_to_structured
nltk.download('punkt')
os.environ['PYTHONHASHSEED'] = '0'  # addition
from datetime import datetime

   

class SCellBOW_pretrain:
    """
    Create the pretrained model from the source dataset.

    Input Arguments:
    - adata_source: the preprocessed scanpy.anndata for source dataset
    - save_dir: name of directory to save the source model
    - vec_size: dimensionality of the embedding vectors. Defaults to 300 for SCellBOW.
    - n_worker: number of worker threads to train the model. For a fully deterministically-reproducible run, limit the model to one worker thread. Defaults to 1 for SCellBOW.
    - iter: Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.
    """
    
    def __init__(self, adata_source, save_dir, vec_size=300, n_worker=1, iter=20):
        self.adata_source = adata_source
        self.corpus_not_shuffled_trn = None
        self.corpus_shuffled_trn = None
        self.save_dir = './{}/'.format(save_dir)
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)
        self.vec_size = vec_size
        self.n_worker = n_worker
        self.iter = iter
        self.run(self.iter, self.vec_size)
        


    # Create corpus
    def wordfreq_docs_list(self, df):
        corpus = []
        row = []
        s = ''
        names = {e: name for e, name in enumerate(df.columns)}
        for i in df.itertuples():
            for e, j in enumerate(i[1:]):
                temp = names[e]
                row += [temp] * int(j)
            corpus.append(row)
            s = ''
            row = []
        
        return corpus

    #Shuffle corpus
    def shuf(self, cnsl):
        corpus_shuffled = []
        random.seed(0)
        for l in tqdm_notebook(range(len(cnsl))):
            random.shuffle(cnsl[l])
            s = ' '.join(cnsl[l])
            corpus_shuffled.append(s)
        print("[",datetime.now(),"]",'Corpus created with size = {}'.format(len(corpus_shuffled)))
        return corpus_shuffled

    #Train Doc2vec
    def doc2vec(self, corpus, iter, model_name, vec_size):
        # tagging docs
        print("[",datetime.now(),"]",'Tagging the corpora.')
        tagged_data = [TaggedDocument(words=word_tokenize(
            _d), tags=[str(i)]) for i, _d in enumerate(corpus)]
        
        print("[",datetime.now(),"]",'All corpuses tagged with length', len(tagged_data))
  
        print("[",datetime.now(),"]",'Inititalize the SCellBOW source model.')
        vector_size=vec_size
        alpha=0.025  
        min_alpha=0.00025
        min_count=1
        window=2
        workers=self.n_worker
        print("[",datetime.now(),"]",'INFO - SCellBOW: vector size =', vec_size)
        print("[",datetime.now(),"]",'INFO - SCellBOW: initial learning rate =', alpha)
        print("[",datetime.now(),"]",'INFO - SCellBOW: min_alpha =', min_alpha)
        print("[",datetime.now(),"]",'INFO - SCellBOW: min_count =', min_count)
        print("[",datetime.now(),"]",'INFO - SCellBOW: number of cpu =', self.n_worker)
    
        model = Doc2Vec(vector_size=vec_size,
                  alpha=0.025,  # initial learning rate
                  min_alpha=0.00025,
                  min_count=1,
                  window=2,
                  workers=self.n_worker,
                  seed=0,
                  dm=1)

        # Building a vocabulary
        print("[",datetime.now(),"]",'Building vocabulary.')
        model.build_vocab(tagged_data, update=False)
        print("[",datetime.now(),"]",'Vocabulary built.')     
        print("[",datetime.now(),"]","Start training the neural network.")
        model.train(tagged_data,
                    total_examples=model.corpus_count,
                    epochs=iter)
        print("[",datetime.now(),"]","Training SCellBOW source model finished.")
        # Save Model
        model.save(self.save_dir + model_name)
        print("[",datetime.now(),"]","Model saved in directory ", self.save_dir)
        return None
   

    def run(self, iter, vec_size):
        print("[",datetime.now(),"]","The path to save directory is" , self.save_dir)
        print("[",datetime.now(),"]","Creating the source model.")
        #rescale the data
        srcdata = self.adata_source.to_df()
        scaler = MinMaxScaler(feature_range=(1, 10))
        scaler.fit(srcdata)
        trainData = scaler.transform(srcdata)
        print("[",datetime.now(),"]","Creating the corpus.")
        trainData = pd.DataFrame(trainData, columns=srcdata.columns)
        # Shuffle the corpus
        self.corpus_not_shuffled_trn = self.wordfreq_docs_list(trainData)
        self.corpus_shuffled_trn = self.shuf(self.corpus_not_shuffled_trn)
        # Train the model
        self.doc2vec(self.corpus_shuffled_trn, iter=self.iter,
                    model_name='source_model', vec_size=vec_size)
        print("[",datetime.now(),"]","Source model created!")
        return None
    
    

    

class SCellBOW_cluster:
    """
    Transfer learning the weights of pre-trained to obtain single-cell embeddings for the target dataset.

    Input Arguments:
    - adata_target: the preprocessed scanpy.anndata for source dataset
    - save_dir: name of directory where the source model is saved
    - resolution: granularity of the leiden clustering. Defaults to 1.0 for SCellBOW.
    - neighbors: number of neighboring data points. Defaults to 15 for SCellBOW.
    - iter: Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.

    
    """
    #Transfer Learning
    def __init__(self,adata_target,save_dir,iter=40,resolution=1.0,neighbors=15):
        self.adata_target = adata_target
        self.adata = None
        self.save_dir = './{}/'.format(save_dir)
        if not os.path.exists(self.save_dir):
             raise Exception("["+str(datetime.now())+"]"+" Source model path not found.")
        self.iter = iter
        self.big_model = Doc2Vec.load(self.save_dir + 'source_model')
        if not os.path.isfile(self.save_dir+'source_model'):
            raise Exception("["+str(datetime.now())+"]"+" Source model not found.")
        self.resolution = resolution
        self.neighbors = neighbors
        
    

    def wordfreq_docs_list(self, df):
        corpus = []
        row = []
        s = ''
        names = {e: name for e, name in enumerate(df.columns)}
        for i in df.itertuples():
            for e, j in enumerate(i[1:]):
                temp = names[e]
                row += [temp] * int(j)
            corpus.append(row)
            s = ''
            row = []
        return corpus

    def shuf(self, cnsl):
        corpus_shuffled = []
        random.seed(0)
        for l in tqdm_notebook(range(len(cnsl))):
            random.shuffle(cnsl[l])
            s = ' '.join(cnsl[l])
            corpus_shuffled.append(s)
        return corpus_shuffled

    def run(self):
        print("[",datetime.now(),"]","The path to save directory is" , self.save_dir)
        print("[",datetime.now(),"]","Begin SCellBOW: transfer learning.")
        self.adata = self.adata_target.copy()
        dstdata = self.adata.to_df()
       
        scaler = MinMaxScaler(feature_range=(1, 10))
        print(scaler.fit(dstdata))
        trainData = scaler.transform(dstdata)
        trainData = pd.DataFrame(trainData, columns=dstdata.columns)

        
        print("[",datetime.now(),"]","Creating the corpus.")
        corpus_not_shuffled = self.wordfreq_docs_list(trainData)
        corpus_shuffled = self.shuf(corpus_not_shuffled)

        
        # tokenizing the corpus
        print("[",datetime.now(),"]",'Tagging the corpora for transfer learning.')
        tagged_data_tl = [TaggedDocument(words=word_tokenize(
            _d), tags=[str(i)]) for i, _d in enumerate(corpus_shuffled)]
        print("[",datetime.now(),"]",'All corpuses tagged with length = {}'.format(len(tagged_data_tl)))

        # Updating the vocabulary
        print("[",datetime.now(),"]",'Updating the vocabulary.')
        self.big_model.build_vocab(tagged_data_tl, progress_per=1000, update=True)
        print("[",datetime.now(),"]",'Vocabulary updated.')

        # Retraining the model
        print("[",datetime.now(),"]","Start transfer learning on the neural network.")
        self.big_model.train(
            tagged_data_tl, total_examples=self.big_model.corpus_count, epochs=self.iter)
        print("[",datetime.now(),"]","Weights of the neural network calibrated.")

        # infer vectors for new corpus
        print("[",datetime.now(),"]","Start infering the vectors for target dataset.")      
        re_vecs_W_new = []
        for c in tqdm_notebook(corpus_shuffled):
            c = c.split(' ')
            self.big_model.random.seed(0)  # addition
            re_vecs_W_new.append(self.big_model.infer_vector(
                c, epochs=self.iter))  # addition
        re_vecs_W_new = np.array(re_vecs_W_new)
        print("[",datetime.now(),"]","Embedding created with shape :", re_vecs_W_new.shape)
              
        self.adata.obsm['SCellBOW_embed'] = re_vecs_W_new
        print("[",datetime.now(),"]","Start leiden clustering at resolution:", self.resolution)  
        
        self.adata.obsm['X_embed'] = sc.tl.pca(
            re_vecs_W_new, svd_solver='arpack')
        sc.pp.neighbors(self.adata, n_neighbors=self.neighbors, use_rep = 'X_embed', random_state=0)
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata, key_added='clusters_' +str(self.resolution), resolution=self.resolution)
        #self.adata.write(self.save_dir+"adata_scellbow.h5ad")
        print("[",datetime.now(),"]","SCellBOW clustering has been successful!")
        return self.adata



class SCellBOW_algebra(SCellBOW_cluster):
    """
    Rank the single cell clusters or subtypes based on their relative aggressiveness.

    Input Arguments:

    - adata_test: the unprocessed scanpy.anndata for single-cell data with the annotation(subtype,cluster) in adata_test.obs
    - adata_train: the anndata for bulk RNAseq gene expression matrix with survival data in adata_train.obs
    - save_dir: name of directory where the source model is saved
    - Type: column from adata_test.obs on which we want to classify (subtype/clusters)
    - bootstrap_samples: number of bootstrap iterations. Defaults to 50 for SCellBOW.
    - split: split on single cell dataset. Defaults to 80:20 split for SCellBOW.
    - unit: type of dataset UMI, TPM, FPKM, etc. Default to UMI for SCellBOW.
    - n_top_features: number of top common highly variables genes in bulk RNAseq and single cell RNAseq datasets. Defaults to 1000 for SCellBOW.
    - iter: Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.
    """
    
    def __init__(self, adata_test, adata_train, save_dir, type, algebra = [], iter=40, bootstrap_samples=50, split=0.2, unit="UMI", n_top_features=1000):
        
        self.adata_test = adata_test  # Load single cell data
        self.adata_train = adata_train  # Load Bulk RNAseq data
        self.save_dir = '{}/'.format(save_dir)
        if not os.path.exists(self.save_dir):
            raise Exception("["+str(datetime.now())+"]"+" Source model path not found.")
        self.iter = iter
        self.big_model = Doc2Vec.load(self.save_dir + 'source_model')
        if not os.path.isfile(self.save_dir+'source_model'):
            raise Exception("["+str(datetime.now())+"]"+" Source model path not found.")
        
        self.Type = type
        if not self.Type in self.adata_test.obs: 
                raise Exception("["+str(datetime.now())+"]"+" Invalid attribute name {}. Please enter a valid attribute!".format(self.Type))
        self.bootstrap_samples = bootstrap_samples
        self.split = split
        self.unit = unit
        self.n_top_features = n_top_features
        self.algebra = algebra
        

    # Create document out of Gene Expression matrix

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
        print("[",datetime.now(),"]","Shape of pseudobulk by individual type:", celltype_avg_vec.shape)
        
        # stack All and type wise Average Vector
        stacked_vec = pd.concat([all_vec, celltype_avg_vec])
        
        # Add combined   
        if self.algebra:
            #check if the list is valid or not
            vector = list(set(self.algebra)) #Take only unique values
            # check if the elements in the list belongs to the dataframe
            for ele in vector:
                if ele not in data['type'].values:
                    raise Exception("["+str(datetime.now())+"]"+" Incorrect entry: This value not exists - "+ str(ele))
            
            if len(vector)<2:
                    raise Exception("["+str(datetime.now())+"]"+"Too few arguments for algebraic operations!")
                    
                    
            # vertically stack the vectors    
            combined_vector = pd.DataFrame()
            for ele in vector:
                combined_vector = pd.concat([combined_vector,pd.DataFrame(celltype_avg_vec.loc[ele]).T])
            # sum of all row
            combined_sum = combined_vector.sum()
            combined_sum.name = '+'.join(vector)
            stacked_vec = pd.concat([stacked_vec,pd.DataFrame(combined_sum).T])
            print("[",datetime.now(),"]","Adding combined vector for algebra:", vector)
       
        # Create scanpy object
        adata_stack = sc.AnnData(stacked_vec)
        adata_stack.obs['type'] = stacked_vec.index.values
        return adata_stack

    def preprocessing(self, adata, unit, top_gene_num):
        adata.var_names_make_unique()  # takes unique genes
        sc.pp.filter_genes(adata, min_cells=20)
        if self.unit == "UMI" or self.unit == "CPM":
            sc.pp.normalize_total(adata, target_sum=1e4)    
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=top_gene_num)
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        return adata

    def run(self):
        print("[",datetime.now(),"]","The path to save directory is" , self.save_dir)
        print("[",datetime.now(),"]","Begin SCellBOW: phenotype algebra.")
        
        
        # Class balance of Single cell based on type
        data = self.class_imbalance(self.adata_test, self.Type)
        
        print("[",datetime.now(),"]","Begin creating pseudobulk on", self.Type)
        
        # Prepare Pseudobulk from single cell
        adata_SCell = self.sc_pseudobulk_anndata(data)
        
        # Combine Bulk and single-cell dataset
        self.adata_train.var_names_make_unique()
        adata_SCell.var_names_make_unique()
        adata_target = self.adata_train.concatenate(
            adata_SCell, join='inner', fill_value=0)
        n_obs, n_vars = adata_target.shape
        if n_vars is None:
            raise Exception("["+str(datetime.now())+"]"+"No common gene between survival and target data!")
        else:
            print("["+str(datetime.now())+"]"+"Common gene between survival and target data={}".format(n_vars))
        
        # Preprocess combined dataset
        adata_target = self.preprocessing(
            adata_target, self.unit, self.n_top_features)
        adata_target.var_names_make_unique()

        # Transfer Learning
        print("[",datetime.now(),"]","Begin transfer learning.")

        # Call ScellBOW_test
        # INherit ScellBOW_clust
        SCellBOW_cluster.__init__(self, adata_target, self.save_dir)
        adata_tl=super().run()  
        
        # Load the embedding
        vector = adata_tl.obsm['SCellBOW_embed']
        print("[",datetime.now(),"]","Finished transfer learning.")

        # length of bulk data
        n = len(self.adata_train.to_df())

        # bulk vector
        bulk_vec = vector[:n]
        bulk_obj = sc.AnnData(bulk_vec)
        bulk_obj.obs = self.adata_train.obs
        #print(bulk_obj.shape)

        # single cell vector
        single_vec = vector[n:]
        single_obj = sc.AnnData(single_vec)
        single_obj.obs = adata_SCell.obs
        #print(single_obj.shape)

        # Survival Analysis
        bulk_data = pd.DataFrame(bulk_vec)
        bulk_data['duration'] = self.adata_train.obs.time.values
        bulk_data['event'] = self.adata_train.obs.status.values
        #print(bulk_data.head(5))

        # Separate Single-cell data for Testing

        single_data = single_obj.to_df()
        single_data['type'] = adata_SCell.obs['type'].values
        #print("Shape:", single_data.shape)
        #print(single_data.head(5))

        # Prepare All celltype vs Rest celltype data
        all_vec = single_data.loc[single_data['type'] == 'All']
        all_vec = all_vec.drop(['type'], axis=1)
        rest_vec = single_data.drop(all_vec.index.values)

        # Survival Analysis
        # to store the predicted score for each patient
        print("[",datetime.now(),"]","Start training the phenotype algebra model.")
        print("[",datetime.now(),"]",'INFO - SCellBOW: Samples in survival data shape =', bulk_data.shape[0])
        print("[",datetime.now(),"]",'INFO - SCellBOW: Pseudobulk samples in target data =', single_data.shape[0])
        print("[",datetime.now(),"]",'INFO - SCellBOW: Descriptor class =', self.Type)
        print("[",datetime.now(),"]",'INFO - SCellBOW: Train:test split = {}:{}'.format(int(((1-self.split)*100)),int((self.split)*100)))
        
        predicted_score = defaultdict(list)
        # Bootstrap on rest of the data
        for j in tqdm_notebook(range(0, self.bootstrap_samples)):  # run boostrapping
            df_train = bulk_data
            print("[",datetime.now(),"]","Epoch {}/{}".format(j,self.bootstrap_samples))

            df_test = df_train.sample(frac=self.split, random_state=j)
            df_train = df_train.drop(df_test.index)

            X_train = df_train.drop(['duration', 'event'], axis=1)
            y_train = convert_to_structured(
                df_train['duration'], df_train['event'])
            #print(X_train.shape, y_train.shape)

            # Train the model
            pr = RandomSurvivalForest(
                random_state=0, n_jobs=-1).fit(X_train, y_train)

            # Score for combined pseudobulk
            X_test = pd.DataFrame(all_vec.values)
            key = "pseudobulk"
            y_predicted = pr.predict(X_test)
            predicted_score[key].append(y_predicted.item(0))

            # Prediction
            for index, vector in rest_vec.iterrows():
                key = "pseudobulk - (" + vector['type']+")"
                test_celltype = pd.DataFrame(vector.drop(['type']).values).T
                sub_test = pd.DataFrame(all_vec.values - test_celltype.values)
                predicted_score[key].append(pr.predict(sub_test).item(0))

        predicted_risk_score = pd.DataFrame(predicted_score)
        print("[",datetime.now(),"]","Risk score prediction complete.")
        print("[",datetime.now(),"]","Calculate median risk score.")          
        median_risk_score = predicted_risk_score.median()
        print("[",datetime.now(),"]","SCellBOW phenotype algebra is complete!")          
        # save the predicted risk score
        return median_risk_score, predicted_risk_score