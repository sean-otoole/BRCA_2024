"""
This script provides utilities for calculating the normalized  ranked AUC for specific gene sets.

It is based on previous work, specificvally the AUCell package written for R.

See R packages here: https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html

Example:
    You can use the functions in this module as follows:
    
    anndata_object = AUCell(anndata_object,signatures_dict,signatures_names gene_lists, AUC_threshold)

Will return a series of normalized gene signature scores that are stored in the anndata_object.
"""

import anndata as ad
import pandas as pd
import numpy as np
from scipy.stats import rankdata


def rank_rows_random(sparse_matrix):
    matrix = sparse_matrix.toarray()

    ranked_matrix = np.zeros_like(matrix, dtype=float)
    for i, row in enumerate(matrix):
        # Get random permutations of the row indices
        random_order = np.random.permutation(len(row))
        # Apply rankdata to the row with the method 'ordinal' to avoid averaging in case of ties
        # Use the random order to break ties
        ranked_row = rankdata(row[random_order], method='ordinal')
        # Undo the random permutation to align with original row order
        ranked_matrix[i, random_order] = ranked_row

    # Rank should be descending, so subtract from the length of the row + 1
    ranked_matrix = len(matrix[0]) + 1 - ranked_matrix
    return ranked_matrix

def calculate_max_auc(current_sig,auc_threshold):   # Calculate max AUC (Ideal case where all genes are in top ranks)
    x_th = np.arange(1, len(current_sig) + 1)   #ideal ranks
    x_th = np.sort(x_th[x_th < auc_threshold])  # ideal ranks less than threshold
    y_th = np.arange(1, len(x_th) + 1)  #cumulative sums for y axis
    max_auc = np.sum(np.diff(np.concatenate(([0], x_th, [auc_threshold]))) * np.concatenate(([0], y_th)))   #integrate across the curve
    return max_auc

def calculate_auc(one_ranking, auc_threshold, max_auc):
    x = np.array(one_ranking)
    # Filter and sort ranks below the threshold
    x = x[x < auc_threshold]
    x = np.sort(x)
    # Create cumulative count
    y = np.arange(1, len(x) + 1)
    # Calculate the actual AUC
    auc = np.sum(np.diff(np.concatenate(([0], x, [auc_threshold]))) * np.concatenate(([0], y)))
    normalized_auc = auc/max_auc
    return normalized_auc

def AUCell(adata,signatures_dict,signatures_names,auc_threshold_percentage = 0.05):
    rank_mat = rank_rows_random(adata.X)
    rank_df = pd.DataFrame(rank_mat, columns = adata.var['gene_ids'].index)   #a ranking dataframe used for all subsequent calculations
    auc_threshold = auc_threshold_percentage*max(rank_df.iloc[0,:])  #lowest rank to consider for top genes
    for sig_name in signatures_names:   #loop through the gene signatures
        current_sig = signatures_dict[sig_name]   # grab the list of genes for each signatures
        current_sig = [gene for gene in current_sig if gene in rank_df.columns]  # filter genes out of signature that are not present in the dataset
        max_auc = calculate_max_auc(current_sig,auc_threshold)
        norm_auc_list = []
        for index, row in rank_df.iterrows():
            sig_rankings = row[current_sig]
            current_auc = calculate_auc(sig_rankings,auc_threshold,max_auc)
            norm_auc_list.append(current_auc)
        adata.obs[sig_name] = norm_auc_list
    return adata
        
        




        
        
