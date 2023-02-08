<center><h1>SCellBOW tutorial </h1></center>
<center> Contributor: Namrata Bhattacharya, Sam Koshy Thomas, Sanket Deshpande</center>
<hr>
<br>
 

# 1. Getting started with SCellBOW
<hr>

## 1.1. Installation Instructions


```python
pip install -i https://test.pypi.org/simple/ SCellBOW0328
```

## 1.2. The base SCellBOW function


```python
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



```python
SCellBOW_cluster(adata_target,save_dir,resolution=1.0,neighbors=15, iter=20,).run()
```

> Transfer learning the weights of pre-trained to obtain single-cell embeddings for the target dataset. 
> #### Input Arguments
> The arguments are as follows:
> - **adata_target:**  the preprocessed scanpy.anndata for target dataset
> - **save_dir:** name of directory where the source model is saved
> - **resolution:** granularity of the leiden clustering. Defaults to 1.0 for SCellBOW. 
> - **neighbors:** number of neighboring data points. Defaults to 15 for SCellBOW. 
> - **iter:** Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.


```python
SCellBOW_algebra(adata_target, adata_surv, save_dir, type,  bootstrap_samples=50, split=0.2, unit="UMI", n_top_features=1000, iter=20).run()
```

> Rank the single cell clusters or subtypes based on their relative aggressiveness.
> #### Input Arguments
> The arguments are as follows:
> - **adata_test:**  the scanpy.anndata for single-cell data with the annotation(subtype,cluster) in *adata_test.obs*
> - **adata_train:**  the anndata for bulk RNAseq gene expression matrix with survival data in *adata_train.obs*
> - **save_dir:** name of directory where the source model is saved
> - **type:** column from *adata_test.obs* on which we want to classify (subtype/clusters).
> - **algebra:** values from column *Type* from *adata_test.obs* which we want to combine (*optional*). 
> - **bootstrap_samples:** number of bootstrap iterations. Defaults to 50 for SCellBOW. 
> - **split:** split on single cell dataset. Defaults to 80:20 split for SCellBOW.
> - **unit:** type of dataset UMI, TPM, FPKM, etc. Default to UMI for SCellBOW. 
> - **n_top_features:** number of top common highly variables genes in bulk RNAseq and single cell RNAseq datasets. Defaults to 1000 for SCellBOW.
> - **iter:** Number of iterations (epochs) over the corpus. Defaults to 20 for SCellBOW.

# Full tutorial 
<hr>

For step-by-step guide on how SCellBOW perform clustering and phenotype algebra on full single-cell dataset.

## 1.  SCellBOW for clustering

* ### 1.1. Import python libraries


```python
import SCellBOW as sb
import scanpy as sc
import matplotlib.pyplot as plt
```

- ###  1.2. Read source input dataset


```python
adata_source = sc.read("/path/to/directory/adata_source.h5ad")
```

- ###  1.3. Preprocessing the count matrix


```python
adata_source.var_names_make_unique()
sc.pp.filter_cells(adata_source, min_genes=200)
sc.pp.filter_genes(adata_source, min_cells=20)

sc.pp.normalize_total(adata_source, target_sum=1e4)
sc.pp.log1p(adata_source)
    
sc.pp.highly_variable_genes(adata_source, n_top_genes = 1000)
adata_source = adata_source[:, adata_source.var.highly_variable]

sc.pp.scale(adata_source, max_value=10)
```

- ###  1.4. Call SCellBOW_pretrain()


```python
sb.SCellBOW_pretrain(adata_source, 'dummy', vec_size=300, n_worker=1, iter=20)
```

- ### 1.5. Read source input dataset


```python
adata_target = sc.read("/path/to/directory/adata_target.h5ad")
```

- ###  1.6. Preprocessing the count matrix


```python
adata_target.var_names_make_unique()
sc.pp.filter_cells(adata_target, min_genes=200)
sc.pp.filter_genes(adata_target, min_cells=20)

sc.pp.normalize_total(adata_target, target_sum=1e4)
sc.pp.log1p(adata_target)
    
sc.pp.highly_variable_genes(adata_target, n_top_genes = 1000)
adata_target = adata_target[:, adata_target.var.highly_variable]

sc.pp.scale(adata_target, max_value=10)
```

- ### 1.7. Call SCellBOW_cluster()


```python
adata = sb.SCellBOW_cluster(adata_target,'dummy').run()
```

- ### 1.8. Visualization


```python
resolution = 1.0
with plt.rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, 
               color='clusters_'+str(resolution), 
               add_outline=True, 
               legend_fontsize=14, 
               legend_fontoutline=2,
               title='UMAP visualisation', 
               size = 50,
               palette=plt.rcParams["axes.prop_cycle"],
              )
```

##  2. SCellBOW for phenotype algebra

- ### 2.1. Import python libraries


```python
import SCellBOW as sb
import scanpy as sc
import matplotlib.pyplot as plt
```

- ###  2.2. Read source input dataset

```python
adata_source = sc.read("/path/to/directory/adata_source.h5ad")
```

- ### 2.3. Read single-cell input dataset


```python
adata_test = sc.read("/path/to/directory/adata_target.h5ad")
```

- ### 2.4. Read bulk survival input dataset


```python
adata_surv = sc.read("/path/to/directory/adata_survival.h5ad")
```
Note that for phenotype algebra the bulk expression dataset has to be in the `AnnData` object format. Ensure that the column names in the `.obs` attribute follow the specific naming convention as below-

|       |time     |status|
| :--:    | :--:      | :--:|
|0_bulk |8.442838 |	1 |
|1_bulk |29.829172|	1 |
|2_bulk |4.533509 |	1 |
|3_bulk |22.996058|	1 |
|4_bulk |28.318003|	0 |

Follow the code snippet to generate the survival data:

```python
import scanpy as sc
import pandas as pd

# Load bulk expression matrix into a Pandas dataframe
df = pd.read_csv('data.csv')

# Load the corresponding survival information into a Pandas dataframe
label = pd.read_csv('label.csv')

# Create an AnnData object from the dataframe
adata = sc.AnnData(df, obs_names=df.index, var_names=df.columns)

# Add sample annotations to the AnnData object
adata.obs = label
adata.obs = adata.obs.rename(columns={'old_colname_time': 'time', 'old_colname_status': 'status'})
```


- ### 2.5. Call SCellBOW_algebra() 


```python
median_score, scores = sb.SCellBOW_algebra(adata_target,
                                           adata_surv,
                                           save_dir ='dummy', 
                                           type='subtype',
                                           algebra = [],
                                           bootstrap_samples=50, 
                                           split=0.2, 
                                           unit="UMI", 
                                           n_top_features=1000, 
                                           iter=20).run()
```

- ### 2.6. Visualization


```python
median_score.sort_values(ascending=True, inplace=True)
scores = scores[median_score.index]
plt.figure(figsize=(4,4))
scores.boxplot(patch_artist=True, notch=True)
plt.xticks(rotation=90, size=10)
plt.show()
```

Please see our **jupyter notebook** examples.
| Notebooks | Description |
| --- | --- |
| [Clusters](RUN-SCellBOW-cluster.ipynb) | Run SCellBOW cluster on pancreas dataset |
| [Phenotype algebra](RUN-SCellBOW-algebra.ipynb) | Run phenotype algebra on molecular subtypes of prostate cancer on a simulated prostate cancer dataset |
| [Phenotype algebra on clusters](RUN-cluster_with_algebra.ipynb) | Run phenotype algebra and SCellBOW clusters on a simulated prostate cancer dataset  |

<br>
Note: Create a folder as the save_dir before running the script.

<br>

## Contributing
Souce code: [Github](https://github.com/cellsemantics/SCellBOW)  
For further information contact debarka@iiitd.ac.in or namrata.bhattacharya@hdr.qut.eu.au. 
