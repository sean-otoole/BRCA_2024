# Start the script and print a message to indicate the process has begun
cat("R script started\n")

# Load necessary libraries
library(Seurat)  # For single-cell RNA-seq analysis
library(dplyr)   # For data manipulation
library(future)  # For parallel computing support

# Set options for parallel processing
options(future.globals.maxSize = 20 * 1024^3)  # Set maximum global memory size to 20 GiB (for large datasets)
plan("multicore", workers = 2)  # Use 2 cores for parallel processing

# Set the plan to use sequential processing (so only one core is used for the SCTransform step)
future::plan("sequential")

# Specify the path to the Seurat object file
objects_path <- paste0(getwd(), '/objects/seurat_objects.rds')

# Load the Seurat object from the specified file path
seurat_objects = readRDS(objects_path)

# Convert 'nFeature_Spatial' to numeric if it's not already in the correct format
seurat_objects[[]]$nFeature_Spatial <- as.numeric(seurat_objects[[]]$nFeature_Spatial)

# Extract the metadata from the Seurat object
metadata <- seurat_objects[[]]

# Create patient-specific dummy variables to account for batch effects during normalization
# The model.matrix function creates a matrix of dummy variables (one for each patient)
sample_dummies <- model.matrix(~ sample - 1, data = metadata)

# Convert the dummy variable matrix to a data frame
sample_dummies_df <- as.data.frame(sample_dummies)

# Add the dummy variables to the Seurat object's metadata
seurat_objects@meta.data <- cbind(seurat_objects@meta.data, sample_dummies_df)

# Specify the variables to regress out during SCTransform (e.g., technical factors and patient-specific batch effects)
current_reg_vars <- c('nCount_Spatial', 'nFeature_Spatial', 'percent.mt')

# Add the patient-specific dummy variables to the list of regressors
current_reg_vars <- append(current_reg_vars, colnames(sample_dummies_df))  # Append patient dummy variable names

# Remove dashes from the column names to avoid potential issues in the regression model
colnames(seurat_objects@meta.data) <- gsub("-", "", colnames(seurat_objects@meta.data))

# Also remove dashes from the names of variables to be regressed out
current_reg_vars <- gsub("-", "", current_reg_vars)

# Perform SCTransform normalization on the Seurat object, regressing out specified variables
# This step normalizes the data, corrects for unwanted technical variation, and selects highly variable genes
seurat_objects <- SCTransform(seurat_objects, 
                        ncells = 3000,                # Use data from 3000 cells for normalization
                        verbose = TRUE,               # Show progress messages
                        vst.flavor = "v2",            # Use the 'v2' variant of the Variance Stabilizing Transformation (VST)
                        conserve.memory = FALSE,      # If set to true will not return corrected umis for all genes
                        return.only.var.genes = FALSE, # return all genes after normalization
                        assay = 'Spatial',            # Specify which assay to use ('Spatial' in this case)
                        vars.to.regress = current_reg_vars,  # Regress out the specified variables
                        clip.range = c(-4, 4),        # Clip the values of the normalized data to this range
                        variable.features.n = 3000)   # Select the top 3000 variable genes for downstream analysis

# Run PCA for dimensionality reduction, using the first 30 principal components (PCs)
seurat_objects <- seurat_objects %>%
  RunPCA(verbose = FALSE, npcs = 30) %>%  # Perform PCA (reduce to 30 PCs)
  
  # Perform UMAP (Uniform Manifold Approximation and Projection) for visualization
  RunUMAP(dims = 1:30, n.neighbors = 20, min.dist = 0.1)  # Use the first 30 PCs for UMAP

# Save the updated Seurat object to the original file path
saveRDS(seurat_objects, objects_path)

# Print a message to indicate the script has completed successfully
cat("R script completed successfully\n")
