import anndata as ad #annotated data matrix, used for generating single-cell objects in python
from helpFuncs import calculate_similarity
import pandas as pd


def cosineSignatures(adata,pb_df):
    pb_df = pb_df
    # get the list of highly variable genes
    highly_variable_genes = adata.var['highly_variable']
    #hvgs as a list
    hvg_list = highly_variable_genes[highly_variable_genes].index.tolist()
    ## not all hvgs are present in both data sets so I am filtering those that are not found in the reference data set
    columns_to_select = hvg_list
    # Check which columns in the list are actually in the DataFrame
    valid_columns = [col for col in columns_to_select if col in pb_df.columns]
    missing_columns = [col for col in columns_to_select if col not in pb_df.columns]
    print("Valid columns:", len(valid_columns))
    print("Missing columns:", len(missing_columns))
    # Only use columns that are in the DataFrame
    valid_columns = [col for col in columns_to_select if col in pb_df.columns]
    # Select the valid columns
    hvg_list = valid_columns
    hvg_adata = adata[:, hvg_list]
    pb_df = pb_df[hvg_list] # need to filter the pseudobulk likst
    #convert to a dense matrix
    hvg_dense = hvg_adata.X.toarray()   
    
    #convert to a dataframe with appropriate column (genes) and row (cellids) names
    df_vis = pd.DataFrame(hvg_dense, index=hvg_adata.obs_names, columns=hvg_adata.var_names)   
    
    for celltype in pb_df.index.tolist():
        current_cos_sim_list = []
        current_sim = celltype + '_cos_sim'
        current_ref = pb_df.loc[celltype]
        for i in range(df_vis.shape[0]):
            current_cell = df_vis.iloc[i]
            current_cos_sim_list.append(calculate_similarity(current_cell,current_ref)) #should probably be an array
        adata.obs[current_sim] = current_cos_sim_list

    return(adata)