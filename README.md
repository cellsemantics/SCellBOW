## SCellBOW : Single-cell RNA-seq analysis and transfer learning using language models

SCellBOW (short for Single Cell Bag of Words) is an unsupervised transfer learning algorithm for clustering scRNA-seq data and phenotype algebra analysis on the clusters using distributed Bag-of-Words. SCellBOW enables the transformation of scRNA-seq data expression matrices to a low-dimensional vector space, where genes are analogous to words, and cells are analogous to documents. Construction of these cell documents enables the direct application of tools from Natural Language Processing (NLP) to the analysis of single-cell sequencing data. SCellBOW shows that NLP in amalgation with transfer learning allows a better representation of cell clusters which can help to decipher the heterogeneity of cells, especially in cancer genomics. SCellBOW facilitates robust identification of malignant clusters and enables stratification of the clusters based on their survival-risk.


![SCellBOW workflow](Data/images/SCellBOW.png)

<!-- \For thorough details, see our paper: [https://www.nature.com/articles/s41467-020-15851-3](https://www.nature.com/articles/s41467-020-15851-3) -->
<br>

## Usage

The [SCellBOW](https://github.com/cellsemantics/SCellBOW) package is an implementation of Single-cell clustering and algebra using language model Doc2vec. With SCellBOW, the user can perform:

- Preprocessing single-cell gene expression profiles.
- Clustering of the target dataset by transfer learning the weights from the pre-trained model trained on the source dataset.
- Visualize cell clustering results and gene expression patterns.
- Perform algebra to stratify the effect of each phenotype on the tumor progression based on relative aggressiveness at an individual phenotype level. 

<br>


# Installation

To install SCellBOW package you must make sure that your python version is 3.8 +. 

### Prerequisites

Python packages:
- pandas = 1.4.3 
- numpy = 1.23.3 
- gensim = 4.2.0  
- scanpy = 1.9.1
- nltk = 3.7
- scikit-learn = 1.1.1  
- scikit-survival = 0.18.0
- imbalanced-learn = 0.9.1
- xgbse = 0.2.3 
- tqdm = 4.64.1
- matplotlib = 3.5.2



### Setup

SCellBOW is included in PyPI, so you can install it by

```bash
pip install SCellBOW
```

If for some reason this doesn't work, you can also download the package from Github and install it locally:

```bash
conda create -n scellbow python=3.8
source activate scellbow
git clone https://github.com/cellsemantics/SCellBOW.git
cd PyPackage/SCellBOW
pip3 install .
```
<br>

### The base SCellBOW functions

```bash
SCellBOW_pretrain(adata_source, save_dir, vec_size=300, n_worker=1, iter=20)
```

> Create the pre-trained model from the source dataset.
> #### Input Arguments
> The arguments are as follows:
> - **adata_source:**  the preprocessed scanpy.anndata for source dataset
> - **save_dir:** name of directory to save the source model
> - **vec_size:** dimensionality of the embedding vectors. Defaults to 300 for SCellBOW. 
> - **n_worker:** number of worker threads to train the model. For a fully deterministically-reproducible run, limit the model to one worker thread. Defaults to 1 for SCellBOW. 
> - **iter:** Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.

```bash
SCellBOW_cluster(adata_target,save_dir,resolution=1.0,neighbors=15, iter=20,).run()
```
> Transfer learning the weights of pre-trained to obtain single-cell embeddings for the target dataset. 
> #### Input Arguments
> The arguments are as follows:
> - **adata_target:**  the preprocessed scanpy.anndata for source dataset
> - **save_dir:** name of directory where the source model is saved
> - **resolution:** granularity of the leiden clustering. Defaults to 1.0 for SCellBOW. 
> - **neighbors:** number of neighboring data points. Defaults to 15 for SCellBOW. 
> - **iter:** Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.

```bash
SCellBOW_algebra(adata_test, adata_train, save_dir, Type='clusters', algebra=["type1","type2"],  bootstrap_samples=50, split=0.2, unit="UMI", n_top_features=1000, iter=20).run()
```
> Rank the single cell clusters or subtypes based on their relative aggressiveness.
> #### Input Arguments
> The arguments are as follows:
> - **adata_test:**  the unprocessed scanpy.anndata for single-cell data with the annotation(subtype,cluster) in *adata_test.obs*
> - **adata_train:**  the anndata for bulk RNAseq gene expression matrix with survival data in *adata_train.obs*
> - **save_dir:** name of directory where the source model is saved
> - **Type:** column from *adata_test.obs* on which we want to classify (subtype/clusters)
> - **algebra:** values from column *Type* from *adata_test.obs* which we want to combine (*optional*). 
> - **bootstrap_samples:** number of bootstrap iterations. Defaults to 50 for SCellBOW. 
> - **split:** split on single cell dataset. Defaults to 80:20 split for SCellBOW.
> - **unit:** type of dataset UMI, TPM, FPKM, etc. Default to UMI for SCellBOW. 
> - **n_top_features:** number of top common highly variables genes in bulk RNAseq and single cell RNAseq datasets. Defaults to 1000 for SCellBOW.
> - **iter:** Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.


<br>

## API example usage

For step-by-step tutorials on how SCellBOW can perform clustering and phenotypic algebra can be found in [Tutorial](tutorial/SCellBOW_tutorial.md).


Here is example usage of SCellBOW in Python:

```bash

import SCellBOW as sb

# List of datasets:
adata_source = [scanpy.AnnData for source dataset]
adata_target = [scanpy.AnnData for source dataset]
adata_bulkseq = [AnnData for bulk RNAseq dataset]


# Creating pre-trained model from source dataset
sb.SCellBOW_pretrain(adata_source, save_dir = 'path/to/model').run()

# Retraining the model with target dataset  
adata_target = sb.SCellBOW_clust(adata_target, save_dir = 'path/to/model').run()

# Predict the risk score for the subtypes in the target dataset
median_score, scores = sb.SCellBOW_algebra(adata_target, adata_bulkseq, save_dir = 'path/to/model').run()

```
<br>

## Contributing
Souce code: [Github](https://github.com/cellsemantics/SCellBOW)  
For further information contact debarka@iiitd.ac.in or namrata.bhattacharya@hdr.qut.eu.au. 











