import scanpy as sc  # Import Scanpy for handling single-cell data
import pandas as pd  # Import pandas for handling data structures and I/O operations
from scipy.io import mmread  # Import mmread to read sparse matrix formats
import decoupler as dc  # Import decoupler for generating pseudobulk data
import os  # For file and directory operations


# Paths to the input files

# Get the current working directory
current_directory = os.getcwd()

count_path = current_directory + '/references/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx' # count matrix
genes_path = current_directory + '/references/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv'  # Gene names
barcodes_path = current_directory + '/references/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv'  # Cell barcodes
metadata_path = current_directory + '/references/Wu_etal_2021_BRCA_scRNASeq/metadata.csv'  # Cell metadata

# Read the sparse count matrix and transpose it (genes as columns, barcodes as rows)
count_matrix = mmread(count_path).tocsr()  # Read the sparse matrix in compressed sparse row format
count_matrix = count_matrix.transpose()  # Transpose the matrix to have genes in columns and barcodes in rows

# Load gene names and barcodes
genes = pd.read_csv(genes_path, header=None, sep='\t')  # Read gene names
barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')  # Read cell barcodes
metadata = pd.read_csv(metadata_path, index_col=0)  # Read metadata (cell type, counts, features)

# Create an AnnData object with the count matrix
adata = sc.AnnData(X=count_matrix)
adata.var['gene_names'] = genes[0].values  # Assign gene names to the AnnData object
adata.obs['barcodes'] = barcodes[0].values  # Assign barcodes to the AnnData object
adata.obs['celltype_major'] = metadata['celltype_major'].values  # Assign major cell type to the AnnData object
adata.obs['nCount_RNA'] = metadata['nCount_RNA'].values  # Assign total RNA count to the AnnData object
adata.obs['nFeature_RNA'] = metadata['nFeature_RNA'].values  # Assign number of RNA features to the AnnData object
adata.obs['celltype_minor'] = metadata['celltype_minor'].values  # Assign minor cell type to the AnnData object
adata.obs['subtype'] = metadata['subtype'].values  # Assign subtype to the AnnData object
adata.var_names = adata.var['gene_names']  # Set gene names as the variable names
adata.obs_names = adata.obs['barcodes']  # Set barcodes as the observation names

# Perform basic preprocessing on the AnnData object:
# - Normalize counts to have a target sum of 10,000
# - Log transform the data
# - Identify highly variable genes (HVGs)
# - Perform PCA for dimensionality reduction
# - Compute neighbors for clustering
# - Compute UMAP embeddings for visualization
sc.pp.normalize_total(adata, inplace=True, target_sum=1e4)  # Normalize the counts
sc.pp.log1p(adata)  # Log-transform the data
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)  # Identify top 2000 highly variable genes
sc.pp.pca(adata)  # Perform PCA for dimensionality reduction
sc.pp.neighbors(adata)  # Compute nearest neighbors for clustering
sc.tl.umap(adata)  # Compute UMAP embeddings for visualization

# Add a column indicating the sample name
adata.obs['sample'] = 'wu'

# Generate pseudobulk data using the decoupler package
# - Group the data by 'celltype_major'
# - Sum the counts for each group
pdata = dc.get_pseudobulk(
    adata,
    sample_col='sample',  # Group by the 'sample' column
    groups_col='celltype_major',  # Group by the 'celltype_major' column
    mode='sum',  # Sum the counts for each group
    min_cells=0,  # Include groups with no cells
    min_counts=0  # Include groups with no counts
)

# Normalize the pseudobulk data and convert it to a DataFrame
sc.pp.normalize_total(pdata, target_sum=1e4, inplace=True)  # Normalize pseudobulk data
pb_df = pd.DataFrame(pdata.X)  # Convert the pseudobulk data to a DataFrame
pb_df.index = pdata.obs['celltype_major']  # Set the row index as the cell type
pb_df.columns = pdata.var['gene_names']  # Set the column names as the gene names

# Save the pseudobulk data as a pickle file for later use
pb_df.to_pickle(os.path.join(os.getcwd(), 'references', 'wu_2021_pseudobulk.pkl'))
