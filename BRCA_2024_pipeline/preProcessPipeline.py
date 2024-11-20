import os
import numpy as np
import pandas as pd
import anndata as ad  # Annotated data matrix, used for generating single-cell objects in Python
import squidpy as sq  # Package for handling spatial transcriptomics datasets
import scanpy as sc  # Toolkit for single-cell analysis
import scvi  # Package for probabilistic modeling in single-cell genomics
import torch  # PyTorch library for deep learning
import celltypist  # Tool for automated cell type annotation
from celltypist import models  # Module for accessing pre-trained models in CellTypist

def processVisium(sample, model):
    ## Import and quality control (QC)
    cwd = os.getcwd()  # Get the current working directory
    current_sample_location = sample
    data_dir = '/tmp/work/Visium/' + current_sample_location  # Path to Visium data
    print(data_dir)  # Debugging: print the data directory
    image_dir = current_sample_location + "/spatial"  # Path to spatial images
    current_sample = current_sample_location.split('/')[1]  # Extract the sample name from the path

    # Read the 10x Visium data
    adata = sq.read.visium(data_dir, library_id=current_sample, source_image_path=image_dir)
    adata.var_names_make_unique()  # Make gene names unique, as some are duplicated (e.g., multiple probes)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")  # Identify mitochondrial genes

    # Calculate quality control metrics (e.g., mitochondrial content)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

    ## Quality filtering
    # Filter out cells with fewer than 600 counts
    sc.pp.filter_cells(adata, min_counts=600)
    # Filter out cells expressing fewer than 500 genes
    sc.pp.filter_cells(adata, min_genes=500)
    # Filter out cells with >10% mitochondrial content
    adata = adata[adata.obs["pct_counts_mt"] < 10].copy()
    print(f"#cells after MT filter: {adata.n_obs}")  # Debugging: print the number of cells after filtering
    # Filter out genes expressed in fewer than 10 locations
    sc.pp.filter_genes(adata, min_cells=10)

    ## Standard processing steps
    sc.pp.normalize_total(adata, inplace=True, target_sum=1e4)  # Normalize gene expression to a total count of 10,000
    sc.pp.log1p(adata)  # Log-transform the data
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=1000)  # Identify the top 1000 highly variable genes
    sc.pp.pca(adata)  # Perform Principal Component Analysis (PCA) for dimensionality reduction
    sc.pp.neighbors(adata)  # Compute the k-nearest neighbors graph
    sc.tl.umap(adata)  # Perform UMAP embedding for visualization
    sc.tl.leiden(adata, key_added="clusters", directed=False, n_iterations=2)  # Cluster cells using the Leiden algorithm

    ## Cell type annotation using CellTypist
    # Predict cell types using the provided model and add majority voting for consensus predictions
    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    adata = predictions.to_adata()  # Update the AnnData object with the predictions

    return adata  # Return the processed AnnData object
