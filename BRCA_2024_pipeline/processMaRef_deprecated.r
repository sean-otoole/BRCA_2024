# Load required libraries
library(Seurat)     # For single-cell RNA-seq analysis and visualization
library(dplyr)      # For data manipulation (e.g., piping and transforming data)
library(gridExtra)  # For arranging multiple plots in a grid
library(ggplot2)    # For creating and customizing plots (ggtitle comes from ggplot2)
library(future)     # For parallel processing to improve computational efficiency

# Set maximum global memory size to 30 GB for future operations (e.g., large data objects)
options(future.globals.maxSize = 50 * 1024^3)  # 100 GiB limit for globals

# Set up parallel processing with 2 workers/cores to speed up computations
#plan("multicore", workers = 1)  # Reduce the number of cores to 2 to reduce memory usage
future::plan("sequential")

# Define the path to the reference dataset (processed single-cell RNA-seq data)
# The reference data is from a public repository and stored locally
ref_data_dir <- paste0(getwd(), '/references/zenodo.11468564/scRNA_data/panB_scRNA_processed_data.rds')

# Load the reference dataset from the specified directory (single-cell RNA-seq data)
ref_data <- readRDS(ref_data_dir)

ref_data@meta.data[is.na(ref_data@meta.data)] <- 0   #get rid of NA values

# subset the data to improve performance
set.seed(123)
subset_cells <- sample(Cells(ref_data), 100000)  # Subset 1000 random cells
ref_data <- subset(ref_data, cells = subset_cells)

# Define a custom color palette with 36 colors for visualizing different cell types
# This will help differentiate the cell types in the UMAP plot
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

# Create the first UMAP plot (before SCTransform) with a title
# This plot visualizes the data before any normalization (SCTransform) has been applied
p1 <- DimPlot(ref_data, group.by = "celltype", label = TRUE, cols = my36colors) + 
  NoLegend() +                          # Remove the legend for clearer visualization
  ggtitle("Before SCTransform")         # Add a title to the plot

# Perform SCTransform normalization, regressing out specified variables (e.g., technical and patient-specific factors)
metadata <- ref_data[[]]

# Create patient-specific dummy variables to account for batch effects
patient_dummies <- model.matrix(~ patient - 1, data = metadata)
patient_dummies_df <- as.data.frame(patient_dummies)

# Add patient dummies to Seurat metadata
ref_data@meta.data <- cbind(ref_data@meta.data, patient_dummies_df)

# Specify variables to regress out (technical factors + patient dummies)
current_reg_vars <- c('percent.mt','DIG.Score1','S.Score','nFeature_RNA','nCount_RNA','G2M.Score')
current_reg_vars <- append(current_reg_vars, colnames(patient_dummies_df))  # Add patient dummies

# Normalize data, regress out unwanted variables, and select top variable genes
ref_data <- SCTransform(ref_data, 
                        ncells = 5000, 
                        verbose = TRUE, 
                        vst.flavor = "v2", 
                        conserve.memory = FALSE,   # If set to true will not return corrected umis for all genes
                        return.only.var.genes = FALSE, 
                        assay = 'RNA',
                        vars.to.regress = current_reg_vars, 
                        clip.range = c(-4, 4), 
                        variable.features.n = 3000) %>% 
  RunPCA(verbose = FALSE, npcs = 30) %>%  # Run PCA for dimensionality reduction
  RunUMAP(dims = 1:30, n.neighbors = 20, min.dist = 0.1)  # Perform UMAP visualization

# Save the transformed reference dataset to an RDS file for future use
# This allows us to re-use the transformed data without having to reprocess it
ma_et_al_path <- paste0(getwd(), '/objects/b_cell_ref_data.rds')
saveRDS(ref_data, ma_et_al_path)

# Create the second UMAP plot (after SCTransform) with a title
# This plot visualizes the data after normalization (SCTransform) has been applied
p2 <- DimPlot(ref_data, group.by = "celltype", label = TRUE, cols = my36colors) + 
  NoLegend() +                          # Remove the legend for clarity
  ggtitle("After SCTransform")          # Add a title to the plot

# Combine the two UMAP plots (before and after SCTransform) into a single figure
# The plots will be arranged side by side (ncol = 2)
combined_plot <- grid.arrange(p1, p2, ncol = 2)

# Define the path to save the combined figure as a PNG file in the 'figures' directory
ma_et_al_path <- paste0(getwd(), '/figures/ma_et_al_umaps.png')

# Save the combined figure as a high-resolution PNG file
# The plot will have dimensions of 14x8 inches with a resolution of 300 DPI
ggsave(ma_et_al_path, plot = combined_plot, width = 14, height = 8, dpi = 300)
