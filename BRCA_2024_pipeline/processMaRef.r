# Load required libraries
library(Seurat)     # For single-cell RNA-seq analysis
library(dplyr)      # For data manipulation (e.g., piping and transforming)
library(gridExtra)  # For arranging multiple plots in a grid
library(ggplot2)    # For creating and customizing plots (ggtitle comes from ggplot2)

# Define the path to the reference dataset (processed single-cell RNA-seq data)
ref_data_dir <- paste0(getwd(), '/references/zenodo.11468564/scRNA_data/panB_scRNA_processed_data.rds')

# Load the reference dataset from the specified directory
ref_data <- readRDS(ref_data_dir)

# Define a custom color palette for the different cell types (36 colors)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

# Create the first UMAP plot (before SCTransform) with a title
# This plot visualizes the data before normalization (SCTransform)
p1 <- DimPlot(ref_data, group.by = "celltype", label = TRUE, cols = my36colors) + 
  NoLegend() +                          # Remove the legend for better clarity
  ggtitle("Before SCTransform")         # Add title to the plot

# Perform SCTransform normalization on the reference data
# This step normalizes the gene expression data and reduces unwanted technical variation
ref_data <- SCTransform(ref_data, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%          # Perform Principal Component Analysis (PCA) on the normalized data
    RunUMAP(dims = 1:30)                 # Run UMAP to perform dimensionality reduction (use the first 30 PCs)

# Save the transformed reference dataset to an RDS file for future use
saveRDS(ref_data, 'b_cell_ref_data.rds')

# Create the second UMAP plot (after SCTransform) with a title
# This plot visualizes the data after normalization (SCTransform)
p2 <- DimPlot(ref_data, group.by = "celltype", label = TRUE, cols = my36colors) + 
  NoLegend() +                          # Remove the legend
  ggtitle("After SCTransform")          # Add title to the plot

# Combine the two UMAP plots into a single figure (side by side) for easy comparison
# The plots will be arranged in two columns (ncol = 2)
combined_plot <- grid.arrange(p1, p2, ncol = 2)

# Save the combined figure as a PNG file with specified dimensions and high resolution
ggsave("ma_et_al_umaps.png", plot = combined_plot, width = 14, height = 8, dpi = 300)

