## SCellBOW : Single-cell RNA-seq analysis and transfer learning using language models

SCellBOW is an unsupervised transfer learning algorithm for clustering scRNA-seq data and phenotype algebra analysis on the clusters using distributed Bag-of-Words. SCellBOW enables the transformation of scRNA-seq data expression matrices to a low-dimensional vector space, where genes are analogous to words, and cells are analogous to documents. Construction of these cell documents enables the direct application of tools from Natural Language Processing (NLP) to the analysis of single-cell sequencing data. SCellBOW shows that NLP in amalgation with transfer learning allows a better representation of cell clusters which can help to decipher the heterogeneity of cells, especially in cancer genomics. SCellBOW facilitates robust identification of malignant clusters and enables stratification of the clusters based on their survival-risk.


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
- pandas = 0.22.0, 
- numpy = 1.16.4, 
- keras = 2.2.4, 
- scipy = 1.0.1, 
- scanpy = 1.3.1+21.g1df151f, 
- anndata = 0.6.20, 
- natsort = 5.2.0, 
- sklearn = 0.19.1



### Setup

SCellBOW is included in PyPI, so you can install it by

```bash
pip install -i https://test.pypi.org/simple/ SCellBOW0328
```

If for some reason this doesn't work, you can also download the package from Github and install it locally:

```bash

```
<br>

## API example usage

For step-by-step tutorials on how SCellBOW can perform clustering and phenotypic algebra can be found in [Tutorial](https://eleozzr.github.io/desc/tutorial.html).


Here is example usage of SCellBOW in Python:

```bash

import SCellBOW as sb

# List of datasets:
adatas = [ list of scanpy.AnnData ]

sb.SCellBOW_pretrain(adata, 
                     save_dir, 
                     vec_size=300, 
                     n_worker=1, 
                     iter=20)


adata_target = sb.SCellBOW_clust(adata_target,
                                 save_dir,iter=40,resolution=1.0,neighbors=15).run()

how to set a 
p=SCellBOW_algebra(adata_test,adata_train,"dd2",Type="clusters",bootstrap_samples=5,n_top_features=50)

t=p.SCellBOW_phenotype_algebra()
```


## Contributing
Souce code: [Github](https://github.com/cellsemantics/SCellBOW)  
For further information contact debarka@iiitd.ac.in or namrata.bhattacharya@hdr.qut.eu.au. 











