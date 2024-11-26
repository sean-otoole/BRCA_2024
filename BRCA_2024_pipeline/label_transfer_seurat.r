# Load required libraries
library(Seurat)     # For single-cell RNA-seq analysis and visualization
library(dplyr)      # For data manipulation (e.g., piping and transforming)
library(future)     # For parallel processing to improve performance

# Set options to control memory size for large objects (20 GiB)
options(future.globals.maxSize = 20 * 1024^3)  # Maximum size limit for globals (set to 130 GB)

# Set up parallel processing plan using 4 CPU cores (adjust based on available cores)
plan("multicore", workers = 4)

# Define the file path for the reference data (processed single-cell RNA-seq data)
ref_data_path <- paste0(getwd(), '/objects/b_cell_ref_data.rds')

# Load the reference dataset from the specified path
ref_data <- readRDS(ref_data_path)

# Define the file path for the Visium dataset (spatial transcriptomics data)
visium_data_path <- paste0(getwd(), '/objects/seurat_objects.rds')

# Load the Visium dataset (spatial transcriptomics data)
visium <- readRDS(visium_data_path)

# Perform PCA (Principal Component Analysis) on the Visium dataset
# This is done to reduce dimensionality and prepare for downstream analysis
visium <- RunPCA(object = visium)

# Split the Visium dataset into a list of Seurat objects based on the "sample" metadata
# Each item in the list will represent a subset of the dataset corresponding to a specific sample
visium_split <- SplitObject(visium, split.by = "sample")

# Loop over each split sample in the Visium dataset to transfer cell type predictions
for (i in seq_along(visium_split)) {
    
    # Extract the current sample from the split dataset
    item <- visium_split[[i]]
    
    # Find transfer anchors between the reference (ref_data) and query (item) Seurat objects
    # Anchors help map the features of the query dataset to the reference dataset
    anchors <- FindTransferAnchors(reference = ref_data, query = item, normalization.method = "SCT", recompute.residuals = FALSE)
    
    # Transfer cell type predictions from the reference to the query using the transfer anchors
    # Use PCA dimensions (1:30) for the transfer process
    predictions.assay <- TransferData(
        anchorset = anchors,                  # Transfer anchors calculated in the previous step
        refdata = ref_data$celltype,          # Reference data cell types
        prediction.assay = TRUE,              # Specify that the result should be stored as a prediction assay
        weight.reduction = item[["pca"]],     # Use PCA reduction for transfer (to speed up computations)
        dims = 1:30                           # Use the first 30 principal components
    )
    
    # Store the predictions in the item (i.e., the current sample Seurat object)
    item[["predictions"]] <- predictions.assay
    
    # Update the list with the modified item containing the predictions
    visium_split[[i]] <- item
}

# Merge the list of Seurat objects (each corresponding to a different sample) back into a single Seurat object
visium <- merge(x = visium_split[[1]], y = visium_split[-1])

# Save the updated Visium dataset with predictions back to disk
saveRDS(visium, visium_data_path)