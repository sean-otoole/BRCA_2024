import re
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import scanpy as sc

def gene_sig_cos_sim_heatmap(condition, sigs, input_dict, order=None, save_path=None):
    """
    Generates a heatmap correlating the cosine similarity scores for specific cell types with a select number of gene signatures.

    Args:
        condition (string): A string used to subset the input dictionary based on sample condition.
        sigs (list): A list of gene signatures to be used in the heatmap.
        input_dict (dict): A dictionary of anndata objects containing data for multiple samples.
        order (list, optional): A list to specify the order of rows and columns within the heatmap. Default is None.
        save_path (str, optional): A string specifying the file path to save the generated plot. Default is None.

    Returns:
        list: List containing the reordered indices of rows and columns from the dendrogram.
    """
    
    # Subset the input dictionary and merge conditions into a single object
    adata_dict = {key: value for key, value in input_dict.items() if condition in key}
    
    # Add a new column 'sample' to each AnnData object to store the sample name
    for sample_name, adata in adata_dict.items():
        adata.obs['sample'] = sample_name  
    
    # Concatenate the selected AnnData objects along the axis=0 (rows), using the sample names as keys
    subset_adata = sc.concat(adata_dict.values(), axis=0, label="sample", keys=adata_dict.keys(), join="outer")
    
    # List of gene signatures of interest
    sigs_of_interest_list = sigs

    # Extract the metadata columns from the combined AnnData object
    meta_columns = subset_adata.obs.columns.tolist()
    
    # Identify columns related to cosine similarity scores
    cos_sim_re = [(i, text, re.search('cos_sim', text)) for i, text in enumerate(meta_columns) if re.search('cos_sim', text)]
    cos_sim_of_interest_list = [sublist[1] for sublist in cos_sim_re]

    # Extract cosine similarity and gene signature columns
    cos_sim_columns = pd.DataFrame({col: subset_adata.obs[col] for col in cos_sim_of_interest_list})
    gene_sig_columns = pd.DataFrame({col: subset_adata.obs[col] for col in sigs_of_interest_list})

    # Ensure both DataFrames have the same number of rows
    if len(cos_sim_columns) != len(gene_sig_columns):
        raise ValueError("The two DataFrames must have the same number of rows.")

    # Fill NaN values with 0 to avoid issues in correlation calculations
    cos_sim_columns.fillna(0, inplace=True)
    gene_sig_columns.fillna(0, inplace=True)

    # Ensure the data types are numeric for correlation calculation
    cos_sim_columns = cos_sim_columns.apply(pd.to_numeric, errors='coerce')
    gene_sig_columns = gene_sig_columns.apply(pd.to_numeric, errors='coerce')

    # Compute pairwise correlations between cosine similarity and gene signature columns
    correlation_result = cos_sim_columns.corrwith(gene_sig_columns)

    # Combine the cosine similarity and gene signature data into a single DataFrame for visualization
    combined_df = pd.concat([cos_sim_columns, gene_sig_columns], axis=1)

    df1 = cos_sim_columns
    df2 = gene_sig_columns

    # Initialize a DataFrame to store the correlation results
    correlation_results = pd.DataFrame(index=df1.columns, columns=df2.columns)

    # Compute pairwise correlations for each pair of columns
    for col1 in df1.columns:
        for col2 in df2.columns:
            correlation_results.loc[col1, col2] = df1[col1].corr(df2[col2])

    # Convert the correlation results to numeric values (in case of non-numeric entries)
    correlation_results = correlation_results.apply(pd.to_numeric)

    # If no custom ordering is provided, generate the default clustermap
    if order is None:
        # Generate the seaborn clustermap
        ax = sb.clustermap(
            correlation_results, annot=None, cmap='coolwarm', center=0,
            linewidths=0, rasterized=True, vmin=-0.3, vmax=0.6, cbar_pos=(0.83, 0.212, 0.03, 0.62)
        )

        # Suppress the row and column dendrograms for cleaner visualization
        ax.ax_row_dendrogram.set_visible(False)
        ax.ax_col_dendrogram.set_visible(False)

        # Remove grid lines inside the heatmap
        ax.ax_heatmap.grid(False)

        # Move the row labels to the left side and ensure they are horizontally aligned
        ax.ax_heatmap.yaxis.set_ticks_position('left')
        ax.ax_heatmap.yaxis.set_tick_params(rotation=0)

        # Set the title for the heatmap
        ax.ax_heatmap.set_title(
            condition + ' cell type cos sim vs Gene signature pairwise correlation',
            pad=20, fontsize=14, loc='center'
        )

        # If a save path is provided, save the plot to the specified location
        if save_path:  
            plt.savefig(save_path, bbox_inches='tight')
            print(f"Plot saved to {save_path}")

        # Return the reordered indices of rows and columns from the dendrogram
        return [ax.dendrogram_row.reordered_ind, ax.dendrogram_col.reordered_ind]

    else:
        # Reorder the correlation results according to the provided order
        correlation_results = correlation_results.iloc[order[0], order[1]]

        # Plot the heatmap using the reordered correlation results
        plt.figure(figsize=(10, 8))
        ax = sb.heatmap(
            correlation_results, annot=None, cmap='coolwarm', center=0,
            linewidths=0, rasterized=True, vmin=-0.3, vmax=0.6
        )

        # Remove grid lines inside the heatmap
        ax.grid(False)

        # Set the title for the heatmap
        plt.title(condition + ' cell type cos sim vs Gene signature pairwise correlation')

        # If a save path is provided, save the plot to the specified location
        if save_path:  
            plt.savefig(save_path, bbox_inches='tight')
            print(f"Plot saved to {save_path}")

    return plt