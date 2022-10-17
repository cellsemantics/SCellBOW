
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
        self.SCellBOW_source(self.iter, self.vec_size)
        

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
   

    def SCellBOW_source(self, iter, vec_size):
        
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