# Script for scRNA-seq reference derived similarity metrics computation

# Load Required Libraries
library(Seurat)      # Single-cell RNA-seq analysis and visualization
library(dplyr)       # Data manipulation (e.g., piping, filtering)
library(gridExtra)   # Arrange multiple plots in a grid
library(ggplot2)     # Data visualization (customized plots)
library(future)      # Parallel processing for computational efficiency
library(transport)   # Earth Mover's Distance (EMD) calculations

# Set Global Options
options(warn = -1)  # Suppress warnings globally (use cautiously)
set.seed(42)        # Set seed for reproducibility
options(future.globals.maxSize = 50 * 1024^3)  # Set max memory for `future` (50 GiB)
future::plan("sequential")  # Use sequential processing to reduce memory usage

# User-Defined Parameters
down_sample <- FALSE  # Toggle down-sampling for testing purposes

#############################################################################################################

# Define Helper Functions
calculate_cosine_similarity <- function(matrix, vector) {
  # Calculate cosine similarity between a vector and columns of a matrix
  vector_normalized <- vector / sqrt(sum(vector^2))  # Normalize vector
  cosine_similarities <- apply(matrix, 2, function(column) {
    column_normalized <- column / sqrt(sum(column^2))  # Normalize column
    sum(column_normalized * vector_normalized)  # Cosine similarity
  })
  return(cosine_similarities)
}

#############################################################################################################

# Load Reference Data
ref_data_dir <- paste0(getwd(), '/references/zenodo.11468564/scRNA_data/panB_scRNA_processed_data.rds')
ref_data <- readRDS(ref_data_dir)

# Down-sample Reference Data (Optional)
if (down_sample) {
  subset_cells <- sample(Cells(ref_data), 50000)  # Randomly sample 50,000 cells
  ref_data <- subset(ref_data, cells = subset_cells)
}

#############################################################################################################

# Split Data into Training and Test Sets
cell_names <- colnames(ref_data)
train_cells <- sample(cell_names, size = floor(0.8 * length(cell_names)), replace = FALSE)
test_cells <- setdiff(cell_names, train_cells)

ref_train <- subset(ref_data, cells = train_cells)
ref_test <- subset(ref_data, cells = test_cells)

cat("Training set cells:", ncol(ref_train), "\n")
cat("Test set cells:", ncol(ref_test), "\n")


#############################################################################################################

# Define Color Palette for Visualization
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
                '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
                '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
                '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', 
                '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')


#############################################################################################################

# Visualize Reference Data (UMAP)
p1 <- DimPlot(ref_data, group.by = "celltype", label = TRUE, cols = my36colors) + 
  ggtitle("Ma et al, 2024 - Reference Dataset")
output_path <- paste0(getwd(), '/figures/', 'ref_umap.png')
ggsave(filename = output_path, plot = p1, width = 8, height = 8, dpi = 300)

#############################################################################################################

# Process and Compare Matrices
ref_train <- FindVariableFeatures(ref_train, selection.method = "dispersion", nfeatures = 7000)
ref_features <- VariableFeatures(ref_train)
visium_path <- paste0(getwd(), '/objects/seurat_objects.rds')
visium_merge <- readRDS(visium_path)

if (down_sample) {
  subset_cells <- sample(Cells(visium_merge), 1000)
  visium_merge <- subset(visium_merge, cells = subset_cells)
}

visium_matrix <- GetAssayData(visium_merge, assay = "SCT", slot = "data")
pseudo_b_ref <- AggregateExpression(ref_train, features = ref_features, assays = "RNA", slot = "data", 
                                     return.seurat = TRUE, group.by = "celltype")
pseudo_bulk_matrix <- GetAssayData(pseudo_b_ref, assay = "RNA", slot = "data")


#############################################################################################################

# Filter for Overlapping Genes
overlapping_genes <- intersect(rownames(visium_matrix), rownames(pseudo_bulk_matrix))
visium_matrix <- visium_matrix[overlapping_genes, ]
pseudo_bulk_matrix <- pseudo_bulk_matrix[overlapping_genes, ]
ref_test_matrix <- GetAssayData(ref_test, assay = "RNA", slot = "data")[overlapping_genes, ]

output_path <- paste0(getwd(), '/results/', 'overlapping_genes.rds')
saveRDS(overlapping_genes, file = output_path)

cat("Number of overlapping genes:", length(overlapping_genes), "\n")


#############################################################################################################

# Calculate Cosine Similarity for Each Cell Type
for (type in colnames(pseudo_bulk_matrix)) {
  visium_merge[[type]] <- calculate_cosine_similarity(visium_matrix, pseudo_bulk_matrix[, type])
}

#############################################################################################################

# Generate a Violin plot for the cos scores within the test set

# Set the celltype factor levels before plotting
dnb_cell_name <- 'B.09.DUSP4+AtM'

# Calculate cosine similarities
dn_b <- pseudo_bulk_matrix[, dnb_cell_name]
dn_b_cell_cos <- calculate_cosine_similarity(ref_test_matrix, dn_b)
ref_test[[dnb_cell_name]] <- dn_b_cell_cos

# Set the celltype as factor and ensure the order
sorting_order <- sort(unique(ref_test$celltype))  # Sorting the unique celltypes
ref_test$celltype <- factor(ref_test$celltype, levels = sorting_order)  # Update factor levels
Idents(ref_test) <- ref_test$celltype  # Synchronize Idents with updated factor levels

# Generate the violin plot
vp <- VlnPlot(ref_test, features = dnb_cell_name, cols = my36colors, pt.size = 0) +
  ggtitle("B.09.DUSP4+AtM Similarity Scores Across Celltypes (Test Set)") +  # Set title
  scale_x_discrete(limits = sorting_order) +  # Ensure x-axis categories are in correct order
  guides(fill = guide_legend(order = 1)) +   # Make legend consistent with x-axis order
  theme(plot.margin = unit(c(1, 1, 1, 2), "cm"))  # Top, Right, Bottom, Left

# Statistical Analysis: Pairwise Wilcoxon Test and Kruskal-Wallis Test
similarity_scores <- FetchData(ref_test, vars = c(dnb_cell_name, "celltype"))

# Separate dnb_cell and other cell types
dnb_cell_scores <- similarity_scores[similarity_scores$celltype == dnb_cell_name, dnb_cell_name]
other_scores <- similarity_scores[similarity_scores$celltype != dnb_cell_name, ]

# Initialize results data frame
results <- data.frame(celltype = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Perform pairwise Wilcoxon rank-sum test
for (cell in unique(other_scores$celltype)) {
  current_scores <- other_scores[other_scores$celltype == cell, dnb_cell_name]
  test_result <- wilcox.test(dnb_cell_scores, current_scores)
  results <- rbind(results, data.frame(celltype = cell, p_value = test_result$p.value))
}

# Adjust for multiple testing using FDR correction
results$adjusted_p_value <- p.adjust(results$p_value, method = "fdr")

# Print results
print(results)

# Perform Kruskal-Wallis test
kruskal_test <- kruskal.test(similarity_scores[[dnb_cell_name]] ~ similarity_scores$celltype)
print(kruskal_test)

# Save results to a file
output_path_results <- paste0(getwd(), '/results/dnb_violin_cell_statistics.csv')
write.csv(results, file = output_path_results, row.names = FALSE)

# Save violin plot
output_path_plot <- paste0(getwd(), '/figures/violin_cos_ref.png')
ggsave(filename = output_path_plot, plot = vp, width = 10, height = 8, dpi = 300)

#############################################################################################################

# Earth Mover's Distance (EMD) Analysis

pre <- subset(visium_merge, subset = sample %in% c('PRE-01','PRE-02','PRE-04','PRE-05'))
post <- subset(visium_merge, subset = sample %in% c('POST-03','POST-05','POST-06','POST-11'))

# Initialize an empty data frame
results_df <- data.frame(
  cell_type = character(),       # For cell type names (character column)
  earth_movers_distance = numeric(),  # For EMD values (numeric column)
  p_value = numeric(),           # For p-values (numeric column)
  stringsAsFactors = FALSE       # Prevent strings from being factors
)

for (type in colnames(pseudo_bulk_matrix)){

pre_data <- pre[[type]]
post_data <- post[[type]]

vector1 <- as.numeric(unlist(pre_data))  # Flatten and convert pre_data
vector2 <- as.numeric(unlist(post_data)) # Flatten and convert post_data

# Compute the observed EMD
emd_obs <- wasserstein1d(vector1, vector2)

# Combine the two datasets
combined_data <- c(vector1, vector2)

# Number of permutations
n_permutations <- 10000
emd_permutations <- numeric(n_permutations)

# Perform permutation test
for (i in 1:n_permutations) {
  # Resample the combined data to create two new groups (same size as original)
  resampled_groups <- sample(combined_data, length(vector1))  # Resample for the "pre" group
  resampled_post <- setdiff(combined_data, resampled_groups)  # The remaining data is the "post" group
  
  # Compute the EMD with resampled data
  emd_permutations[i] <- wasserstein1d(resampled_groups, resampled_post)
}

# Calculate p-value (proportion of permuted EMDs greater than observed EMD)
p_val <- mean(emd_permutations >= emd_obs)

results_df <- rbind(results_df, data.frame(
cell_type = type,
earth_movers_distance = emd_obs,
p_value = p_val))
}

# Save Results and Visualize

output_path <- paste0(getwd(), '/results/', 'eds_results_df.rds')
saveRDS(results_df, file = output_path)

emd_plot <- ggplot(results_df, aes(x = cell_type, y = earth_movers_distance, color = cell_type)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) +
  labs(title = "Earth Mover's Distance Comparison", 
       x = "Cell Type", 
       y = "EMD") +
  scale_color_manual(values = my36colors) +  # Apply custom color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(-0.1, "inches")
  )

output_path <- paste0(getwd(), '/figures/', 'emd_plot.png')
ggsave(filename = output_path, plot = emd_plot, width = 8, height = 8, dpi = 300)

#############################################################################################################

# Save Updated Visium Object
saveRDS(visium_merge, file = visium_path)
