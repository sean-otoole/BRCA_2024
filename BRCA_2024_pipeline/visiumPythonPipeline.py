# === Import Necessary Libraries ===
# Standard Libraries
import os  # For handling file paths and directory operations
import sys  # For system-specific parameters and functions
from datetime import datetime  # For handling dates and times
import time  # For time measurements
import subprocess  # For running shell commands
import warnings  # To handle warnings
import pickle  # For saving and loading Python objects

# Data Analysis and Visualization Libraries
import pandas as pd  # For handling tabular data
import numpy as np  # For numerical operations
from scipy.stats import zscore  # For Z-score normalization
from sklearn.decomposition import PCA  # For dimensionality reduction (PCA)
from sklearn.cluster import KMeans  # For K-Means clustering
import seaborn as sb  # For advanced plotting
import matplotlib.pyplot as plt  # For basic plotting

# Specialized Libraries
import pyreadr  # For reading R data files (e.g., .Rds, .RData)
import squidpy as sq  # For spatial transcriptomics analysis
import scanpy as sc  # For single-cell data analysis
from celltypist import models  # For cell type annotation

# Custom Scripts
from applyImageMask import applyMask  # Custom function for image masking
from applyCosineScores import cosineScores  # Custom function for cosine similarity
from getRSigs import getRSigs  # Custom function to retrieve gene signatures
from pyAUCell import AUCell  # Custom function for AUCell scoring
from geneSigCosSimHeatmap import gene_sig_cos_sim_heatmap  # Custom function for generating a specific correaltional heatmap
from preProcessPipeline import processVisium

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# === Define Input Data and Variables ===
# List of Visium sample identifiers
samples = [
    '2018838/PRE-01', '2018838/PRE-02', '2018838/POST-03', '2018838/POST-05',
    '2018839/PRE-04', '2018839/PRE-05', '2018839/POST-06', '2018839/POST-11'
]

# Download the CellTypist model
print("Downloading CellTypist model...")
models.download_models(force_update=True, model='Immune_All_High.pkl')
model = models.Model.load(model='Immune_All_High.pkl')  # Load the model for use

# Create a dictionary to store results for each sample
sample_names = [os.path.basename(sample) for sample in samples]
visium_samples = {name: None for name in sample_names}

# === Step 1: Generate Pseudobulk Reference ===
print("Generating pseudobulk reference...")
start_time = time.time()
subprocess.run(['python', 'acquirePseudoBulkRef.py'])  # Run the pseudobulk script
elapsed_time = time.time() - start_time
print(f"Pseudobulk reference generated in {elapsed_time:.2f} seconds.")

# Load the pseudobulk reference
pb_df = pd.read_pickle(os.path.join(os.getcwd(), 'references', 'wu_2021_pseudobulk.pkl'))

# === Step 2: Retrieve Gene Signatures ===
signatures_dict, signature_groups_list = getRSigs()

# === Step 3: Process Samples ===
start_time = datetime.now()

for i, sample in enumerate(samples):
    print(f"\nProcessing sample: {sample_names[i]}")

    # Process Visium data
    visium_samples[sample_names[i]] = processVisium(sample, model)
    print(f"Cosine similarities applied for {sample_names[i]}...")
    visium_samples[sample_names[i]] = cosineScores(visium_samples[sample_names[i]], pb_df)
    
    print(f"Pathology mask applied for {sample_names[i]}...")
    visium_samples[sample_names[i]] = applyMask(visium_samples[sample_names[i]], sample_names[i])
    
    print(f"Gene signature scores applied for {sample_names[i]}...")
    visium_samples[sample_names[i]] = AUCell(visium_samples[sample_names[i]], signatures_dict, signature_groups_list)

end_time = datetime.now()
execution_time = end_time - start_time
print(f"Processing completed in {execution_time}.")

# Save the processed samples for future use
with open('visium_samples.pkl', 'wb') as file:
    pickle.dump(visium_samples, file)

# === Step 4: Compare Bulk RNA-Seq and Visium Samples ===
# Initialize PCA
pca = PCA(n_components=2)

# Aggregate signature scores for each sample
signature_scores_df = pd.DataFrame(columns=signature_groups_list, index=sample_names)
for sample in sample_names:
    current_data = visium_samples[sample]
    for signature in signature_groups_list:
        signature_scores_df.at[sample, signature] = np.mean(current_data.obs[signature])

# Load bulk gene signatures from an RDS file
current_directory = os.getcwd()
bulk_sigs_path = os.path.join(current_directory, 'gene_signatures', 'all_signature_scores.Rds')
bulk_sigs = pyreadr.read_r(bulk_sigs_path)[None]  # Read and extract the DataFrame

# Normalize signature scores using Z-score normalization
signature_scores_df = signature_scores_df.apply(pd.to_numeric, errors='coerce')
signature_scores_df_z = signature_scores_df.apply(zscore)

# Combine bulk and Visium scores
combined_sigs = pd.concat([bulk_sigs, signature_scores_df_z])

# Apply PCA to the combined dataset
pca_result = pca.fit_transform(combined_sigs)
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'], index=combined_sigs.index)

# Perform K-Means clustering
kmeans = KMeans(n_clusters=2, random_state=123)
pca_df['cluster'] = kmeans.fit_predict(pca_df[['PC1', 'PC2']])

# === Step 5: Visualize Results ===
# Set up color palette
palette = sb.color_palette("hsv", 2)

# Create scatter plot
plt.figure(figsize=(10, 6))
scatter_plot = sb.scatterplot(
    x='PC1', y='PC2', hue='cluster', palette=palette, data=pca_df,
    s=100, edgecolor='w', alpha=0.7
)

# Highlight specific samples
highlight_samples = ['POST-06', 'PRE-01']
for sample in highlight_samples:
    idx = pca_df.index.get_loc(sample)
    scatter_plot.scatter(
        pca_df.PC1.iloc[idx], pca_df.PC2.iloc[idx],
        color='red', s=100, edgecolor='w', alpha=0.7
    )

# Add labels with a 45-degree angle
for i in range(pca_df.shape[0]):
    scatter_plot.text(
        pca_df.PC1.iloc[i], pca_df.PC2.iloc[i], 
        pca_df.index[i], fontsize=7, ha='right', rotation=-45
    )

# Finalize and save plot
plt.title('PCA Plot with K-Means Clustering')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Cluster')
plt.grid(False)

# Save the plot
figures_dir = os.path.join(current_directory, 'figures')
os.makedirs(figures_dir, exist_ok=True)
plt.savefig(os.path.join(figures_dir, 'pca_kmeans_clustering_highlighted_rotated_labels.png'))

# Define all possible cell type categories for consistent ordering
cell_type_categories = ['Epithelial cells', 'Endothelial cells', 'Fibroblasts', 
                         'B cells', 'Plasma cells', 'Macrophages', 'pDC']

# Initialize an empty DataFrame to store the relative frequencies of cell types for each sample
combined_df = pd.DataFrame(index=cell_type_categories)

# Loop through each sample dataset
for sample in sample_names:
    # Extract the dataset for the current sample
    current_data_set = visium_samples[sample]
    
    # Get the 'majority_voting' column, which contains the assigned cell type for each observation
    # current_counts = current_data_set.obs['majority_voting']
    current_counts = pd.Series(current_data_set.obs['majority_voting'].astype('category'))
    
    # Create a frequency table for the cell types, ensuring all categories are included
    frequency_table = current_counts.value_counts().reindex(cell_type_categories, fill_value=0)
    
    # Normalize the frequency table to obtain relative frequencies
    relative_frequencies = frequency_table / frequency_table.sum()
    
    # Add the relative frequencies for this sample to the combined DataFrame
    combined_df[sample] = relative_frequencies

# Print the resulting DataFrame to verify its structure
print(combined_df)

# Create a simplified summary by combining certain categories
simple_cateogories = {
    'epi': combined_df.iloc[0:1, :].sum(axis=0),       # Sum frequencies of 'Epithelial cells'
    'all_others': combined_df.iloc[1:, :].sum(axis=0) # Sum frequencies of all other cell types
}

# Convert the simplified categories into a new DataFrame
simple_df = pd.DataFrame(data=simple_cateogories)

# Transpose the DataFrame so that categories are rows and samples are columns
simple_df = simple_df.transpose()

# Rearrange the columns so that pre and post samples are grouped with eachother
simple_df = simple_df[['PRE-01', 'PRE-02', 'PRE-04', 'PRE-05', 
                       'POST-03', 'POST-05', 'POST-06', 'POST-11']]

# Plot a stacked bar graph for the proportions of cell types across samples
ax = combined_df.T.plot(kind='bar', stacked=True, figsize=(12, 8), colormap='Set1')  # Higher contrast colormap

# Remove grid lines for cleaner visualization
plt.grid(False)

# Add a title and axis labels to the plot
plt.title('Proportions of Different Cell Types Across Samples')
plt.xlabel('Samples')
plt.ylabel('Proportion')

# Add a legend, placing it outside the plot to avoid overlap
plt.legend(title='Cell Types', bbox_to_anchor=(1.02, 1), loc='upper left')

# Adjust the layout to ensure elements fit well and prevent clipping
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leaves space for the legend on the right

# Save the figure to a specified directory with minimal whitespace
plt.savefig(os.path.join(figures_dir, 'cell_type_dist_celltypist.png'))

#correlate the cosine similarities of specific cell types with the gene signatures
pre_corr_plot_location = os.getcwd() + '/figures/' + 'pre_correlations.png'
post_corr_plot_location = os.getcwd() + '/figures/' + 'post_correlations.png'

pre_order = gene_sig_cos_sim_heatmap('PRE',signature_groups_list,visium_samples,save_path=pre_corr_plot_location)
gene_sig_cos_sim_heatmap('POST',signature_groups_list,visium_samples,pre_order,save_path=post_corr_plot_location)

# save all the anndata objects for metadata transfer to seurat objects within R
for sample in visium_samples.keys():
    current_object = visium_samples[sample]
    current_path = os.getcwd() + "/objects/" + sample + ".h5ad"
    current_object.write(current_path)

