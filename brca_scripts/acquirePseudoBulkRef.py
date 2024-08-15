import scanpy as sc
import pandas as pd
from scipy.io import mmread
import decoupler as dc

count_path = '/tmp/work/Visium/sean/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx'
genes_path = '/tmp/work/Visium/sean/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv'
barcodes_path = '/tmp/work/Visium/sean/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv'
metadata_path = '/tmp/work/Visium/sean/Wu_etal_2021_BRCA_scRNASeq/metadata.csv'

count_matrix = mmread(count_path).tocsr()
count_matrix = count_matrix.transpose()  #ensures that genes are in columns and rows are barcodes

genes = pd.read_csv(genes_path, header=None, sep='\t')
barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')
metadata = pd.read_csv(metadata_path, index_col=0)

adata = sc.AnnData(X=count_matrix)
adata.var['gene_names'] = genes[0].values  # Assuming the gene names are in the first column
adata.obs['barcodes'] = barcodes[0].values  # Assuming the barcodes are in the first column
adata.obs['celltype_major'] = metadata['celltype_major'].values  
adata.obs['nCount_RNA'] = metadata['nCount_RNA'].values 
adata.obs['nFeature_RNA'] = metadata['nFeature_RNA'].values
adata.obs['celltype_minor'] = metadata['celltype_minor'].values 
adata.obs['subtype'] = metadata['subtype'].values 
adata.var_names = adata.var['gene_names']
adata.obs_names = adata.obs['barcodes']

# normalization, hvg, embeddings and umap 
sc.pp.normalize_total(adata, inplace=True, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

adata.obs['sample'] = 'wu'

pdata = dc.get_pseudobulk(
    adata,
    sample_col = 'sample',
    groups_col='celltype_major',
    #layer='counts',
    mode='sum',
    min_cells=0,
    min_counts=0
)

sc.pp.normalize_total(pdata, target_sum=1e4, inplace=True)  #normalize prior to converting to a dataframe
pb_df = pd.DataFrame(pdata.X)
pb_df.index = pdata.obs['celltype_major']
pb_df.columns = pdata.var['gene_names']
pb_df.to_pickle('wu_2021_pseudobulk.pkl')


