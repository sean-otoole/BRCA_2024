import anndata as ad  # Annotated data matrix, used for generating single-cell objects in Python
from helpFuncs import calculate_similarity  # Import the function to calculate cosine similarity
import pandas as pd  # Library for handling dataframes

def cosineScores(adata, pb_df):
    pb_df = pb_df  # Pseudobulk DataFrame containing reference data
    
    # Get the list of highly variable genes (HVGs)
    highly_variable_genes = adata.var['highly_variable']
    # Create a list of the highly variable genes that are marked as True
    hvg_list = highly_variable_genes[highly_variable_genes].index.tolist()

    # Not all HVGs are present in both datasets, so we filter those that are not found in the reference dataset (pb_df)
    columns_to_select = hvg_list
    # Check which columns in the list are actually present in the pseudobulk DataFrame
    valid_columns = [col for col in columns_to_select if col in pb_df.columns]
    missing_columns = [col for col in columns_to_select if col not in pb_df.columns]
    print("Valid columns:", len(valid_columns))  # Debugging: print the number of valid columns
    print("Missing columns:", len(missing_columns))  # Debugging: print the number of missing columns

    # Only keep columns that are present in both the data and the pseudobulk DataFrame
    hvg_list = valid_columns
    
    # Subset the AnnData object (adata) to only include the highly variable genes present in the pseudobulk data
    hvg_adata = adata[:, hvg_list]
    # Filter the pseudobulk DataFrame to only include the valid genes
    pb_df = pb_df[hvg_list]  # Filter pseudobulk data based on the selected HVGs
    
    # Convert the sparse matrix to a dense array for further analysis
    hvg_dense = hvg_adata.X.toarray()   

    # Convert the dense matrix to a DataFrame with the appropriate gene names as columns and cell IDs as rows
    df_vis = pd.DataFrame(hvg_dense, index=hvg_adata.obs_names, columns=hvg_adata.var_names)   
    
    # For each cell type in the pseudobulk DataFrame
    for celltype in pb_df.index.tolist():
        current_cos_sim_list = []  # Initialize a list to store the cosine similarity values
        current_sim = celltype + '_cos_sim'  # Create a label for the current cell type's cosine similarity
        current_ref = pb_df.loc[celltype]  # Get the reference pseudobulk data for the current cell type
        
        # Compare each cell's gene expression to the reference pseudobulk profile
        for i in range(df_vis.shape[0]):
            current_cell = df_vis.iloc[i]  # Select the gene expression profile of the current cell
            current_cos_sim_list.append(calculate_similarity(current_cell, current_ref))  # Calculate cosine similarity
            
        # Add the calculated cosine similarities as a new column in the AnnData object
        adata.obs[current_sim] = current_cos_sim_list

    return adata  # Return the updated AnnData object with cosine similarity scores
