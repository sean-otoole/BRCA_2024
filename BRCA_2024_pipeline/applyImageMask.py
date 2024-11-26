import SimpleITK as sitk  # For reading and processing medical images
import matplotlib.pyplot as plt  # For visualization
import numpy as np  # For numerical operations
import squidpy as sq  # For handling spatial transcriptomics data sets
import scanpy as sc  # For single-cell RNA-seq data analysis
import os  # For file and directory operations

def applyMask(adata, sample_name):
    """
    Applies a mask to spatial transcriptomics data to filter out unwanted spots 
    and generates visualizations before and after applying the mask.

    Parameters:
    - adata: AnnData object containing spatial transcriptomics data.
    - sample_name: Name of the sample, used to locate the corresponding mask file.

    Returns:
    - filtered_adata: AnnData object filtered by the mask.
    """
    # Get the current working directory
    current_directory = os.getcwd()

    # Construct the path to the mask image
    mask_path = current_directory + '/image_masks/' + sample_name + '_Mask.png'

    # Read the mask image and convert it to a NumPy array
    image = sitk.ReadImage(mask_path)
    image_array = sitk.GetArrayFromImage(image)

    # Binarize the image: set pixel values > 0 to 1, otherwise 0
    binarized_image = np.where(image_array > 0, 1, 0)

    # Scale the spatial coordinates using the sample-specific scaling factor
    coordinates = adata.obsm['spatial'] * adata.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']

    # Initialize a list to track which spots overlap with the mask
    overlap_vector = []
    for coord in coordinates:
        x, y = coord

        # Skip invalid coordinates (NaN values)
        if np.isnan(x) or np.isnan(y):
            overlap_vector.append(False)
            continue

        # Round the coordinates to integer values
        x_rounded, y_rounded = int(round(x)), int(round(y))

        # Check if the rounded coordinates are within the mask bounds
        if (0 <= x_rounded < binarized_image.shape[1]) and (0 <= y_rounded < binarized_image.shape[0]):
            # Append True if the corresponding mask pixel is 1 (inside the mask)
            overlap_vector.append(binarized_image[y_rounded, x_rounded] == 1)
        else:
            # Append False if the coordinates are outside the mask bounds
            overlap_vector.append(False)

    # Convert the overlap vector to a NumPy array
    overlap_vector = np.array(overlap_vector)

    # Filter the spatial coordinates and the AnnData object based on the mask
    filtered_coordinates = coordinates[overlap_vector]
    filtered_adata = adata[overlap_vector, :]

    # Extract original spatial coordinates for visualization
    spatial_coords = adata.obsm['spatial'] * adata.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']

    # Create the figures directory if it doesn't exist
    output_dir = "figures"
    os.makedirs(output_dir, exist_ok=True)

    # Define the file path for saving the figure
    file_path = os.path.join(output_dir, f"{sample_name}_segmentation_plots.png")

    # Create a 2x2 grid for subplots
    fig, axs = plt.subplots(2, 2, figsize=(16, 16))

    # Subplot 1: Display the mask image
    axs[0, 0].imshow(image_array, cmap='gray')
    axs[0, 0].set_title('Mask Image')
    axs[0, 0].axis('off')  # Hide axes for better clarity

    # Subplot 2: Plot the cell type distribution after applying the mask
    sq.pl.spatial_scatter(filtered_adata, color=["majority_voting"], ax=axs[0, 1])

    # Subplot 3: Plot the original spatial coordinates (pre-segmentation)
    axs[1, 0].scatter(spatial_coords[:, 0], spatial_coords[:, 1], s=10, c='blue', alpha=0.6)
    axs[1, 0].invert_yaxis()  # Invert the y-axis to match the image orientation
    axs[1, 0].set_title(f'Pre-segmentation ({sample_name})')
    axs[1, 0].set_xlabel('X Coordinate')
    axs[1, 0].set_ylabel('Y Coordinate')

    # Capture the axis limits of the pre-segmentation plot
    x_limits = axs[1, 0].get_xlim()
    y_limits = axs[1, 0].get_ylim()

    # Subplot 4: Plot the filtered spatial coordinates (post-segmentation)
    spatial_coords_filtered = (
        filtered_adata.obsm['spatial'] * adata.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']
    )
    axs[1, 1].scatter(spatial_coords_filtered[:, 0], spatial_coords_filtered[:, 1], s=10, c='blue', alpha=0.6)
    axs[1, 1].invert_yaxis()
    axs[1, 1].set_title(f'Post-segmentation ({sample_name})')
    axs[1, 1].set_xlabel('X Coordinate')
    axs[1, 1].set_ylabel('Y Coordinate')

    # Apply the captured axis limits to the post-segmentation plot
    axs[1, 1].set_xlim(x_limits)
    axs[1, 1].set_ylim(y_limits)

    # # Adding labels to each subplot outside the axes in the top-left corner
    # labels = ['a', 'b', 'c', 'd']  # Subpanel labels
    # for ax, label in zip(axs.flat, labels):
    #     ax.text(1.05, 1.05, label, transform=ax.transAxes, 
    #             fontsize=16, fontweight='bold', va='top', ha='left')

    # Adding labels outside the axes but in the top-left corner
    labels = ['a', 'b', 'c', 'd']  # Subpanel labels
    for ax, label in zip(axs.flat, labels):
        ax.text(-0.1, 1.05, label, transform=ax.transAxes, 
                fontsize=16, fontweight='bold', va='top', ha='left')
        
    # Adjust spacing to accommodate the labels outside the axes
    plt.subplots_adjust(left=0.1, right=0.85, top=0.85, bottom=0.1)

    # Adjust layout to avoid overlap between subplots
    plt.tight_layout()

    # Save the figure to the specified file path
    plt.savefig(file_path, dpi=300)  # Save with 300 DPI for high-quality output

    # Close the plot to free up memory
    plt.close(fig)

    # Return the filtered AnnData object
    return filtered_adata