"""
This script provides a function for providing a mask vector to which will isolate the tumour specific region within the slide.

Will generate a pdf showing the results of the tumour segmentation. The pdf is stored in the data folder.

Additionally, it will apply the tumour segmentation to the adata file and return a subsetted adata object.
"""

import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
import os


# Load the RGB image
def retrieveTumourMask(image_path,adata):
    ### image path should be for the hires h&e image
    
    sample_name = image_path.split('/')[5]
    
    rgb_image = sitk.ReadImage(image_path)
    
    # Convert RGB to grayscale
    if rgb_image.GetNumberOfComponentsPerPixel() == 3:
        gray_image = sitk.VectorIndexSelectionCast(rgb_image, 0, sitk.sitkFloat32)
        gray_image = sitk.RescaleIntensity(gray_image, outputMinimum=0, outputMaximum=255)
    else:
        gray_image = sitk.RescaleIntensity(rgb_image, outputMinimum=0, outputMaximum=255)
    
    # Apply Gaussian blur
    smoothing_filter = sitk.SmoothingRecursiveGaussianImageFilter()
    smoothing_filter.SetSigma(7.0)
    blurred_image = smoothing_filter.Execute(gray_image)
    
    # Apply Otsu Threshold
    threshold_filter = sitk.OtsuThresholdImageFilter()
    threshold_filter.SetInsideValue(1)  # Value for foreground pixels
    threshold_filter.SetOutsideValue(0) # Value for background pixels
    thresholded_image = threshold_filter.Execute(blurred_image)
    
    # Apply Binary Opening to remove small objects
    opening_filter = sitk.BinaryOpeningByReconstructionImageFilter()
    opening_filter.SetKernelRadius(100)  # Adjust radius to control the size of particles to remove
    opened_image = opening_filter.Execute(thresholded_image)
    
    # Apply Binary Closing to fill small holes
    closing_filter = sitk.BinaryClosingByReconstructionImageFilter()
    closing_filter.SetKernelRadius(50)  # Adjust radius to control the size of holes to fill
    closed_image = closing_filter.Execute(opened_image)
    
    # Convert images to numpy arrays for plotting
    rgb_image_array = sitk.GetArrayFromImage(rgb_image)
    gray_image_array = sitk.GetArrayFromImage(gray_image)
    blurred_image_array = sitk.GetArrayFromImage(blurred_image)
    thresholded_image_array = sitk.GetArrayFromImage(thresholded_image)
    opened_image_array = sitk.GetArrayFromImage(opened_image)
    closed_image_array = sitk.GetArrayFromImage(closed_image)
    
    # Check the shape of the RGB image array
    print("Shape of RGB image array:", rgb_image_array.shape)
    
    # Reorder RGB image array to have channels as the last dimension if necessary
    if rgb_image_array.ndim == 3 and rgb_image_array.shape[1] == 3:
        rgb_image_array = np.transpose(rgb_image_array, (1, 2, 0))
    
    # Ensure the RGB array shape is correct
    print("Corrected shape of RGB image array:", rgb_image_array.shape)
    
    # Plot the images
    fig, axs = plt.subplots(1, 6, figsize=(36, 6))
    
    # Show the original RGB image
    axs[0].imshow(rgb_image_array)
    axs[0].set_title('Original RGB Image')
    axs[0].axis('off')
    
    # Show the grayscale image
    axs[1].imshow(gray_image_array, cmap='gray')
    axs[1].set_title('Grayscale Image')
    axs[1].axis('off')
    
    # Show the blurred image
    axs[2].imshow(blurred_image_array, cmap='gray')
    axs[2].set_title('Blurred Image')
    axs[2].axis('off')
    
    # Show the thresholded image
    axs[3].imshow(thresholded_image_array, cmap='gray')
    axs[3].set_title('Thresholded Image')
    axs[3].axis('off')
    
    # Show the opened image
    axs[4].imshow(opened_image_array, cmap='gray')
    axs[4].set_title('Opened Image')
    axs[4].axis('off')
    
    # Show the closed image
    axs[5].imshow(closed_image_array, cmap='gray')
    axs[5].set_title('Closed Image')
    axs[5].axis('off')

    dir_name = "/tmp/work/Visium/BRCA_2024/figures" 
    pdf_path = os.path.join(dir_name, sample_name + '_tumour_segmentation' + '.' + 'pdf')
    
    # Save the figure as a PDF
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close(fig)

    # Assuming coordinates are in the following format:
    coordinates = adata.obsm['spatial'] * adata.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']
    binarized_image = closed_image_array
    
    overlap_vector = []
    
    for coord in coordinates:
        x, y = coord
        
        # Check if either x or y is NaN
        if np.isnan(x) or np.isnan(y):
            overlap_vector.append(False)
            continue
        
        # Round the coordinates to the nearest integer
        x_rounded = round(x)
        y_rounded = round(y)
        
        # Ensure the rounded coordinates are within the image bounds
        if 0 <= x_rounded < binarized_image.shape[1] and 0 <= y_rounded < binarized_image.shape[0]:
            if binarized_image[y_rounded, x_rounded] == 1:
                overlap_vector.append(True)
            else:
                overlap_vector.append(False)
        else:
            # Handle out-of-bounds coordinates (append False)
            overlap_vector.append(False)
    
    
    filtered_coordinates = coordinates[overlap_vector]

    filtered_adata = adata[overlap_vector, :]

    return(filtered_adata)

    