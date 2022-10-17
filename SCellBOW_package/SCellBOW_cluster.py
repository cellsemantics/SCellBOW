import time
import pickle
from tqdm import tqdm
import matplotlib.pyplot as plt
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from sklearn.preprocessing import MinMaxScaler
from nltk.tokenize import word_tokenize
import numpy as np
import pandas as pd
import scanpy as sc
import random
import nltk
import os
nltk.download('punkt')

os.environ['PYTHONHASHSEED'] = '0'  # addition


class SCellBOW_clust:
    #Transfer Learning
    def __init__(self):
        pass
 

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

    def SCellBOW_target(self,adata_target, save_dir, iter=40, resolution=1.0, neighbors=15):
        self.adata_target = adata_target
        self.save_dir = "./"+save_dir+"/"
        if not os.path.exists(self.save_dir):
            raise Exception("Train Model not Found")
        self.iter = iter
        self.big_model = Doc2Vec.load(self.save_dir + 'source_model')
        if not os.path.isfile(self.save_dir+'source_model'):
            raise Exception("Source Model not found")
        self.resolution = resolution
        self.neighbors = neighbors
        
        
        
        
        
        print("SCellBOW_target")
        dstdata = self.adata_target.to_df()
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
        self.adata_target.obsm['SCellBOW_embed'] = re_vecs_W_new
        self.adata_target.obsm['X_embed'] = sc.tl.pca(
            re_vecs_W_new, svd_solver='arpack')

        # sc.pp.neighbors(self.adata_target, n_neighbors=self.neighbors, random_state=0)
        # print(self.adata_target)
        # sc.tl.umap(self.adata_target)

        # sc.tl.leiden(self.adata_target, key_added='clusters_' +str(self.resolution), resolution=self.resolution)
        # print(self.adata_target)

        self.adata_target.write(self.save_dir+"adata_scellbow.h5ad")
        # with plt.rc_context({'figure.figsize': (5, 5)}):
            # sc.pl.umap(self.adata_target, color='clusters_'+str(self.resolution),show=True)

        return self.adata_target
