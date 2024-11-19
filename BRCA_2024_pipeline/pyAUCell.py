"""
This script provides utilities for calculating the normalized ranked AUC for specific gene sets.

It is based on previous work, specifically the AUCell package written for R.

See R packages here: https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html

Example:
    You can use the functions in this module as follows:
    
    anndata_object = AUCell(anndata_object, signatures_dict, signatures_names, gene_lists, AUC_threshold)

Will return a series of normalized gene signature scores that are stored in the anndata_object.
"""

import anndata as ad
import pandas as pd
import numpy as np
from scipy.stats import rankdata

def rank_rows_random(sparse_matrix):
    """
    Convert a sparse matrix to a dense matrix and calculate ranked values row by row.
    Rankings are randomized in the case of ties to avoid bias.

    Args:
        sparse_matrix: A sparse matrix representation of gene expression data.

    Returns:
        A matrix with ranked values for each row, ranked in descending order.
    """
    matrix = sparse_matrix.toarray()
    ranked_matrix = np.zeros_like(matrix, dtype=float)

    for i, row in enumerate(matrix):
        # Randomize the row indices to break ties
        random_order = np.random.permutation(len(row))
        # Apply ranking, using 'ordinal' to avoid averaging ties
        ranked_row = rankdata(row[random_order], method='ordinal')
        # Revert the random permutation to restore the original order
        ranked_matrix[i, random_order] = ranked_row

    # Convert rankings to descending order
    ranked_matrix = len(matrix[0]) + 1 - ranked_matrix
    return ranked_matrix

def calculate_max_auc(current_sig, auc_threshold):
    """
    Calculate the maximum possible AUC for a given gene signature.

    Args:
        current_sig: List of genes in the signature.
        auc_threshold: The rank threshold to consider for the top genes.

    Returns:
        The maximum possible AUC value.
    """
    x_th = np.arange(1, len(current_sig) + 1)  # Ideal ranks for the genes
    x_th = np.sort(x_th[x_th < auc_threshold])  # Filter ranks below the threshold
    y_th = np.arange(1, len(x_th) + 1)  # Cumulative sums for the y-axis
    # Calculate the area under the curve (AUC) by integration
    max_auc = np.sum(np.diff(np.concatenate(([0], x_th, [auc_threshold]))) * np.concatenate(([0], y_th)))
    return max_auc

def calculate_auc(one_ranking, auc_threshold, max_auc):
    """
    Calculate the normalized AUC for a specific ranking of genes.

    Args:
        one_ranking: The ranks of the genes in the current cell/sample.
        auc_threshold: The rank threshold for considering top genes.
        max_auc: The maximum possible AUC for normalization.

    Returns:
        The normalized AUC value.
    """
    x = np.array(one_ranking)
    x = x[x < auc_threshold]  # Filter ranks below the threshold
    x = np.sort(x)  # Sort the filtered ranks
    y = np.arange(1, len(x) + 1)  # Cumulative count for y-axis
    # Calculate actual AUC by integration
    auc = np.sum(np.diff(np.concatenate(([0], x, [auc_threshold]))) * np.concatenate(([0], y)))
    normalized_auc = auc / max_auc  # Normalize the AUC
    return normalized_auc

def AUCell(adata, signatures_dict, signatures_names, auc_threshold_percentage=0.05):
    """
    Compute AUCell scores for a set of gene signatures and store them in an AnnData object.

    Args:
        adata: AnnData object containing gene expression data.
        signatures_dict: Dictionary mapping signature names to lists of genes.
        signatures_names: List of signature names to calculate AUC for.
        auc_threshold_percentage: Percentage threshold for determining top-ranked genes.

    Returns:
        AnnData object with normalized AUC scores added to `adata.obs`.
    """
    # Rank all genes for each cell/sample using random tie-breaking
    rank_mat = rank_rows_random(adata.X)
    # Create a DataFrame of ranked data, with genes as columns
    rank_df = pd.DataFrame(rank_mat, columns=adata.var['gene_ids'].index)
    # Calculate the threshold rank based on the percentage
    auc_threshold = auc_threshold_percentage * max(rank_df.iloc[0, :])

    # Loop through each gene signature to calculate AUCell scores
    for sig_name in signatures_names:
        print('\n')
        print('Current sig:', sig_name)
        # Retrieve the list of genes for the current signature
        current_sig = signatures_dict[sig_name]
        print('Original length is', len(current_sig))
        # Filter out genes not present in the dataset
        current_sig = [gene for gene in current_sig if gene in rank_df.columns]
        print('After filtering the length is', len(current_sig))
        print('\n')

        # Calculate the maximum AUC for normalization
        max_auc = calculate_max_auc(current_sig, auc_threshold)
        norm_auc_list = []

        # Compute the normalized AUC for each cell/sample
        for index, row in rank_df.iterrows():
            sig_rankings = row[current_sig]  # Get ranks for genes in the signature
            current_auc = calculate_auc(sig_rankings, auc_threshold, max_auc)
            norm_auc_list.append(current_auc)

        # Store the normalized AUC scores in the AnnData object
        adata.obs[sig_name] = norm_auc_list

    return adata
