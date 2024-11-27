library(Seurat)
options(warn = -1)
library(sceasy)
library(reticulate)

dataset_paths = c('2018838/PRE-01','2018838/PRE-02','2018838/POST-03','2018838/POST-05',
           '2018839/PRE-04','2018839/PRE-05','2018839/POST-06','2018839/POST-11')

dataset_paths <- lapply(dataset_paths, function(x){paste('/tmp/work/Visium/', x, sep = '')})

seurat_objects <- lapply(dataset_paths, function(x) {Load10X_Spatial(unlist(x))})   #creates a list of seurat objects
sample_names <- lapply(dataset_paths, function(x){tail(unlist(strsplit(unlist(x), split = "/")),n=1)})   #gets the sample name for metadata assignment

# commented this out since I've already filtered with python script preprocess
# seurat_objects <- lapply(seurat_objects, function(x){subset(x,subset = nCount_Spatial > 2000)})

apply_sample_names <- function(x,y){
    x[[]]$sample = y
    return(x)
}
seurat_objects <- mapply(apply_sample_names, seurat_objects, sample_names)  #add sample name information to each seurat object
seurat_objects <- merge(x = seurat_objects[[1]], y = seurat_objects[-1])  #merges object across layers
seurat_objects <- PercentageFeatureSet(seurat_objects, pattern="^MT-", col.name="percent.mt")

# Initialize an empty data frame to store results
combined_metadata <- data.frame()
current_number <- 0

# Loop through unique samples
for (sample in unique(seurat_objects[[]]$sample)) {
  current_number <- current_number + 1
  current_path <- paste0(getwd(), '/objects/', sample, '.h5ad')
  output_file_path <- paste0(getwd(), '/objects/', sample, '.rds')
  
  # Convert the file from h5ad to Seurat format
  sceasy::convertFormat(current_path, from = "anndata", to = "seurat", outFile = output_file_path)
  
  # Read the converted Seurat object
  current_seurat_object <- readRDS(output_file_path)
  
  # Extract metadata
  current_meta <- current_seurat_object[[]]
  
  # Modify rownames by appending "_<current_number>"
  rowname_modifier <- paste0("_", current_number)
  rownames(current_meta) <- paste0(rownames(current_meta), rowname_modifier)
  
  # Concatenate vertically
  combined_metadata <- rbind(combined_metadata, current_meta)
}

subset_seurat_object <- subset(seurat_objects, cells = rownames(combined_metadata))

seurat_objects <- subset(seurat_objects, cells = rownames(combined_metadata))

# Check if the row names match
if (identical(rownames(seurat_objects[[]]), rownames(combined_metadata))) {
    print("Meta-dataframes row names match, proceeding with combining metadata")
    
    # Combine metadata
    combined_df <- cbind(seurat_objects[[]], combined_metadata)
    
    # Replace the metadata in the Seurat object
    seurat_objects[[]] <- combined_df
} else {
    print("Row names do not match. Check the order of cells.")
    # Optional: inspect row names
    print(head(rownames(subset_seurat_object[[]])))
    print(head(rownames(combined_metadata)))
}

object_path <- paste0(getwd(),'/objects/','seurat_objects.rds')

saveRDS(seurat_objects,object_path)

# Run the sctransform script
sc_transform_path <- paste0(getwd(), '/sctransform.r')

# Execute the R script using system2 and capture the output and error messages
output <- system2("Rscript", args = sc_transform_path, wait = TRUE, stdout = TRUE, stderr = TRUE)

# Capture the exit status of the executed R script
exit_status <- attr(output, "status")

# Treat NULL exit status as non-critical and proceed
if (is.null(exit_status) || exit_status == 0) {
    # If exit status is NULL or 0, consider the script executed successfully
    print("sctransform.r script executed successfully or exit status is NULL (non-critical).")
} else {
    # If the exit status is non-zero, it indicates failure
    print(paste("sctransform.r script execution failed with exit status:", exit_status))
    cat("Error Message:\n", output)
    
    # Optionally, write the error log to a file for further inspection
    error_log_path <- paste0(getwd(), '/error_log_sctransform.txt')
    write(output, file = error_log_path)
    cat("The error log has been saved to:", error_log_path, "\n")
}

# Check for potential issues with stderr (if any errors were returned)
if (length(output) > 0) {
    cat("Script Output:\n", output, "\n")
} else {
    cat("No output from the script.\n")
}

# Read the Seurat objects
seurat_objects <- readRDS(object_path)

# Run the script to generate the reference data set (Ma et al. 2024)
# Path to the script
process_Ma_Ref_path <- paste0(getwd(), '/processMaRef.r')

# Execute the R script using system2 and capture the output and error messages
output <- system2("Rscript", args = process_Ma_Ref_path, wait = TRUE, stdout = TRUE, stderr = TRUE)

# Capture the exit status of the executed R script
exit_status <- attr(output, "status")

# Treat NULL exit status as non-critical and proceed
if (is.null(exit_status) || exit_status == 0) {
    # If exit status is NULL or 0, consider the script executed successfully
    print("ProcessMaRef.r script executed successfully or exit status is NULL (non-critical).")
} else {
    # If the exit status is non-zero, it indicates failure
    print(paste("ProcessMaRef.r script execution failed with exit status:", exit_status))
    cat("Error Message:\n", output)
    
    # Optionally, write the error log to a file for further inspection
    error_log_path <- paste0(getwd(), '/error_log_processMaRef.txt')
    write(output, file = error_log_path)
    cat("The error log has been saved to:", error_log_path, "\n")
}

# Check for potential issues with stderr (if any errors were returned)
if (length(output) > 0) {
    cat("Script Output:\n", output, "\n")
} else {
    cat("No output from the script.\n")
}

# Path to the R script you want to execute
label_transfer_path <- paste0(getwd(), '/label_transfer_seurat.r')

# Execute the R script using system2 and capture the output and error messages
output <- system2("Rscript", args = label_transfer_path, wait = TRUE, stdout = TRUE, stderr = TRUE)

# Capture the exit status of the executed R script
exit_status <- attr(output, "status")

# Treat NULL exit status as non-critical and proceed
if (is.null(exit_status) || exit_status == 0) {
    # If exit status is NULL or 0, consider the script executed successfully
    print("label_transfer_seurat.r script executed successfully or exit status is NULL (non-critical).")
} else {
    # If the exit status is non-zero, it indicates failure
    print(paste("label_transfer_seurat.r script execution failed with exit status:", exit_status))
    cat("Error Message:\n", output)
    
    # Optionally, write the error log to a file for further inspection
    error_log_path <- paste0(getwd(), '/error_log_label_transfer.txt')
    write(output, file = error_log_path)
    cat("The error log has been saved to:", error_log_path, "\n")
}

# Check for potential issues with stderr (if any errors were returned)
if (length(output) > 0) {
    cat("Script Output:\n", output, "\n")
} else {
    cat("No output from the script.\n")
}

# Calculate and add the cosine similarity scores
process_cos_ma_ref <- paste0(getwd(), '/cos_sim.r')

# Execute the R script using system2 and capture the output and error messages
output <- system2("Rscript", args = process_cos_ma_ref, wait = TRUE, stdout = TRUE, stderr = TRUE)

# Capture the exit status
exit_status <- attr(output, "status") # Retrieve the exit status from the system2 output attributes

# Treat NULL exit status as non-critical and proceed
if (is.null(exit_status) || exit_status == 0) {
    # If exit status is NULL or 0, consider the script executed successfully
    print("cos_sim.r script executed successfully or exit status is NULL (non-critical).")
} else {
    # If the exit status is non-zero, it indicates failure
    print(paste("cos_sim.r script execution failed with exit status:", exit_status))
    cat("Error Message:\n", output)
    
    # Optionally, write the error log to a file for further inspection
    error_log_path <- paste0(getwd(), '/error_log_cos_sim.txt')
    write(output, file = error_log_path)
    cat("The error log has been saved to:", error_log_path, "\n")
}

# Check for potential issues with stderr (if any errors were returned)
if (length(output) > 0) {
    cat("Script Output:\n", output, "\n")
} else {
    cat("No output from the script.\n")
}

