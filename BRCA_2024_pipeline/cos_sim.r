library(Seurat)
library(dplyr)
library(ggplot2)

calculate_cosine_similarity <- function(matrix, vector) {
  
  # Normalize the vector
  vector_norm <- sqrt(sum(vector^2))
  vector_normalized <- vector / vector_norm
  
  # Initialize a vector to store cosine similarities for each column
  cosine_similarities <- numeric(ncol(matrix))
  
  # Loop through each column in the matrix and calculate cosine similarity
  for (i in 1:ncol(matrix)) {
    
    # Get the column
    column <- matrix[, i]
    
    # Normalize the column
    column_norm <- sqrt(sum(column^2))
    column_normalized <- column / column_norm
    
    # Calculate cosine similarity as dot product of normalized vectors
    cosine_similarities[i] <- sum(column_normalized * vector_normalized)
  }
  
  return(cosine_similarities)
}

ref_path <- paste0(getwd(), '/objects/b_cell_ref_data.rds')

data = readRDS(ref_path)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

# Create the plot
dim_plot <- DimPlot(data, group.by = "celltype", label = TRUE, cols = my36colors) + NoLegend()

fig_path <- paste0(getwd(), '/figures/ma_reference_dimplot.png')

# Save as PNG
ggsave(fig_path, plot = dim_plot, width = 8, height = 6, dpi = 300)

## subset the b cell data set
unique_celltypes <- unique(data[[]]$celltype)
b_cell_types_only <- unique_celltypes[c(1:9,12,13)]
b_ref <- subset(data, subset = celltype %in% b_cell_types_only)

## select the highly variable features
b_ref <- FindVariableFeatures(b_ref, selection.method = "vst", nfeatures = 2000)
b_features <- VariableFeatures(b_ref)

## generate a psuedobulk dataframe where the rows are the HVGs and the columns are cell type
pseudo_b_ref <- AggregateExpression(b_ref, features = b_features, assays = "RNA", return.seurat = T, group.by = "celltype")
b_ref_matrix <- GetAssayData(pseudo_b_ref, assay = "RNA", slot = "data")

## import the seurat processed data
visium_path <- paste0(getwd(), '/objects/seurat_objects.rds')

visium_merge <- readRDS(visium_path)
visium_matrix <- GetAssayData(visium_merge, assay = "SCT", slot = "scale.data")
overlapping_genes <- intersect(rownames(visium_matrix),rownames(b_ref_matrix))
print(length(overlapping_genes))  # will need to expand the number genes for the SCT transformation

## subset the matrices
visium_matrix_overlap <- visium_matrix[overlapping_genes,]
b_ref_matrix_overlap <- b_ref_matrix[overlapping_genes,]

dn_b <- b_ref_matrix_overlap[,'B.09.DUSP4+AtM'] 
naive_b <- b_ref_matrix_overlap[,'B.01.TCL1A+naiveB'] 
just_b <- b_ref_matrix_overlap[,'B.02.IFIT3+B'] 


## calculate cosine similarity for each location and add it to the seurat objects
dn_b_cell_cos <- calculate_cosine_similarity(visium_matrix_overlap,dn_b)
naive_b_cell_cos <- calculate_cosine_similarity(visium_matrix_overlap,naive_b)
just_b_cell_cos <- calculate_cosine_similarity(visium_matrix_overlap,just_b)

visium_merge$dn_b_cos <- dn_b_cell_cos
visium_merge$naive_b_cell_cos <- naive_b_cell_cos
visium_merge$just_b_cell_cos <- just_b_cell_cos

saveRDS(visium_merge,visium_path)