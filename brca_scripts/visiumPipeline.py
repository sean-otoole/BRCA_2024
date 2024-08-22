import os
import numpy as np
import pandas as pd
import anndata as ad  #annotated data matrix, used for generating single-cell objects in python 
import squidpy as sq #package for handling spatial transcriptomics data sets
import scanpy as sc
import scvi
import torch #pytorch
import celltypist
from celltypist import models

def processVisium(sample, model):
    ##import and qc
    cwd = os.getcwd()
    current_sample_location = sample
    data_dir = cwd.strip('BRCA_2024')+current_sample_location
    print(data_dir)
    image_dir = current_sample_location+"/spatial"
    current_sample = current_sample_location.split('/')[1]
    adata = sq.read.visium(data_dir, library_id=current_sample,source_image_path=image_dir)  #reads the 10x visium data
    adata.var_names_make_unique()   # several genes are duplicated (multiple probes?) this removes the duplicates
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

    #quality filtering
    sc.pp.filter_cells(adata, min_counts=600) ## pp = preprocess, filter coordinates with less than 600 reads.
    sc.pp.filter_cells(adata, min_genes=500)  ## filter with less than 500 genes expressed.
    adata = adata[adata.obs["pct_counts_mt"] < 10].copy()
    #print(f"#cells after MT filter: {adata.n_obs}")
    sc.pp.filter_genes(adata, min_cells=10)  ## filter genes expressed in less than 10 locations.

    #standard processing steps
    sc.pp.normalize_total(adata, inplace=True, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=1000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters", directed=False, n_iterations=2)

    # getting models for cell typist
    # in the future this should be placed elsewhere and only run once

    predictions = celltypist.annotate(adata, model = model, majority_voting = True)
    adata = predictions.to_adata()    
    return adata