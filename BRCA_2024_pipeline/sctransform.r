cat("R script started\n")
library(Seurat)
library(dplyr) 
library(future)
options(future.globals.maxSize = 20 * 1024^3)  #15 GiB
plan("multicore", workers = 2)

future::plan("sequential")

objects_path <- paste0(getwd(), '/objects/seurat_objects.rds')
seurat_objects = readRDS(objects_path)

# Perform SCTransform normalization, regressing out specified variables (e.g., technical and patient-specific factors)
metadata <- seurat_objects[[]]

# Create patient-specific dummy variables to account for batch effects
sample_dummies <- model.matrix(~ sample - 1, data = metadata)
sample_dummies_df <- as.data.frame(sample_dummies)

# Add patient dummies to Seurat metadata
seurat_objects@meta.data <- cbind(seurat_objects@meta.data, sample_dummies_df)

# Specify variables to regress out (technical factors + patient dummies)
current_reg_vars <- c('nCount_Spatial','nFeature_Spatial','percent.mt')
current_reg_vars <- append(current_reg_vars, colnames(sample_dummies_df))  # Add patient dummies

# Normalize data, regress out unwanted variables, and select top variable genes
ref_data <- SCTransform(seurat_objects, 
                        ncells = 3000, 
                        verbose = TRUE, 
                        vst.flavor = "v2", 
                        conserve.memory = TRUE, 
                        return.only.var.genes = TRUE, 
                        assay = 'RNA',
                        vars.to.regress = current_reg_vars, 
                        clip.range = c(-4, 4), 
                        variable.features.n = 3000) %>% 
  RunPCA(verbose = FALSE, npcs = 30) %>%  # Run PCA for dimensionality reduction
  RunUMAP(dims = 1:30, n.neighbors = 20, min.dist = 0.1)  # Perform UMAP visualization

saveRDS(seurat_objects, objects_path)
cat("R script completed successfully\n")