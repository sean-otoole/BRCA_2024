options(warn = -1)  # Suppress warnings globally if needed
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(dplyr)
})

spatial_metrics <- function(current_metric, gene = FALSE){
    
    object_path <- paste0(getwd(),'/objects/','seurat_objects.rds')
    seurat_objects <- readRDS(object_path)
    
    # List of Seurat objects, you can replace these with your actual Seurat objects
    samples <- c('PRE-01', 'PRE-02', 'PRE-04', 'PRE-05',
                 'POST-03', 'POST-05', 'POST-06', 'POST-11')
    
    # Create an empty list to store the plots
    plot_list <- list()
    
    # Fetch feature data from the first Seurat object to determine min and max values for the color scale
    first_seurat_obj <- subset(seurat_objects, subset = sample %in% c(samples[1]))
    feature_data <- FetchData(first_seurat_obj, vars = current_metric)
    min_cut <- round(min(feature_data[[current_metric]], na.rm = TRUE),1)
    max_cut <- round(max(feature_data[[current_metric]], na.rm = TRUE),1)
    mid_point <- round(min_cut + (max_cut - min_cut)/2,1)
    
    # Loop through each sample and plot the spatial data
    for (i in 1:length(samples)) {
        sample_name <- samples[i]
        # Subset Seurat object based on sample name (ensure sample is a metadata column)
        seurat_obj <- subset(seurat_objects, subset = sample %in% c(sample_name))
        #Create SpatialFeaturePlot with custom color scale range
        P <- SpatialFeaturePlot(seurat_obj, features = current_metric, pt.size = 2, alpha = 0.8)
        p <- P + ggplot2::scale_fill_gradient(low = "blue", high = "yellow", limits = c(min_cut, max_cut), breaks = c(min_cut, mid_point, max_cut),guide = guide_colorbar(
            label = TRUE,        # keep the labels
            title = NULL,        # Optional: remove the title of the color bar
            barwidth = 10,       # Adjust the width of the color bar
            barheight = 1        # Adjust the height of the color bar
          )) + 
        ggtitle(paste(sample_name))  +  # Add sample name as title
        theme(plot.title = element_text(hjust = 0.5))  + # Center the title# Add sample name as title
        theme(plot.margin=unit(c(0,0.1,0,0.1), "cm")) +
        theme(
        legend.position = "top",  # Moves the legend to the top of the plot
        legend.title = element_text(size = 10),  # Adjust legend title size
        legend.text = element_text(size = 8),  # Adjust legend text size
        legend.key.height = unit(0.5, "lines"))  # Adjust the height of the legend keys
        #Append the plot to the list
        plot_list[[i]] <- p
    }
    
    # Arrange the plots into a grid layout (2 rows, 4 columns) without default title
    grid_plot_upper <- grid.arrange(
        grobs = plot_list, 
        ncol = 4,
        top = NULL  # Remove the default title, we add it manually
    )
    
    
    pre <- subset(seurat_objects, subset = sample %in% c('PRE-01','PRE-02','PRE-04','PRE-05'))
    post <- subset(seurat_objects, subset = sample %in% c('POST-03','POST-05','POST-06','POST-11'))
    
    column_name <- current_metric
    
    options(repr.plot.width = 7, repr.plot.height = 7)  # This will make the plot twice as wide
    
    # Combine data into a dataframe, accounting for different sizes

    if (gene) {
        # if a gene is being, use the expression data
        pre_data_gene_exp <- GetAssayData(pre,layer = c("scale.data"))
        post_data_gene_exp <- GetAssayData(post,layer = c("scale.data"))
        pre_data <- pre_data_gene_exp[column_name,]
        post_data <- post_data_gene_exp[column_name,]
    } else {
            # otherwise, use the metadata
        pre_data <- pre[[column_name]]
        post_data <- post[[column_name]]
    }
    
    df <- data.frame(value = rbind(pre_data, post_data),
                    group = factor(c(rep("Pre", length(unlist(pre_data))),rep("Post", length(unlist(post_data))))))
    
    colnames(df)[1] <- "value"
    
    # # Plot with ggplot
    dist_plot <- ggplot(df, aes(x = value, fill = group)) + 
    geom_density(alpha = 0.5) + 
    labs(title = "Density Plot of Pre and Post Groups", x = "Value", y = "Density") +
    scale_fill_manual(values = c("red", "blue")) +
    theme_minimal()
    
    # Define the feature to analyze
    current_obs <- current_metric
    
    # Define the post and pre sample groups
    post_samples <- c('POST-03', 'POST-05', 'POST-06', 'POST-11')
    pre_samples <- c('PRE-01', 'PRE-02', 'PRE-04', 'PRE-05')
    
    # Subset Seurat object based on sample group and extract the feature data for "post" and "pre" samples
    post_values <- sapply(post_samples, function(current_sample) {
      # Subset based on the sample metadata and fetch the feature data
      subset_data <- subset(seurat_objects, subset = sample == current_sample)
      mean(FetchData(subset_data, vars = current_obs)[[current_obs]], na.rm = TRUE)
    })
    
    pre_values <- sapply(pre_samples, function(current_sample) {
      subset_data <- subset(seurat_objects, subset = sample == current_sample)
      mean(FetchData(subset_data, vars = current_obs)[[current_obs]], na.rm = TRUE)
    })
    
    # Perform a t-test
    t_test_result <- t.test(pre_values, post_values, var.equal = FALSE)
    cat("T-statistic:", round(t_test_result$statistic, 4), "p-value:", round(t_test_result$p.value, 4), "\n")
    
    # Prepare data for plotting
    data <- data.frame(
      Sample = rep(c('Pre', 'Post'), times = c(length(pre_values), length(post_values))),
      Value = c(pre_values, post_values)
    )
    
    # Create a plot with strip plot, horizontal mean line, no grid, ticks inside the plot, and uniform font size
    means_plot <- ggplot(data, aes(x = Sample, y = Value)) +
      geom_jitter(width = 0.1, alpha = 0.7, color = 'skyblue') +  # Blue circles for individual samples
      stat_summary(fun = "mean", geom = "crossbar", color = "red", size = 1.5, width = 0.2, fatten = 1) + # Red horizontal line for means
      labs(
        title = paste('Pre vs Post', current_obs, 'Averages\nT-test p-value:', round(t_test_result$p.value, 4)),
        x = 'Sample Group',
        y = paste('Average', current_obs)
      ) +
      theme_minimal() +  # Minimal theme (removes some default elements)
      theme(
        panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(size = 1, color = "black"),  # Add axis lines
        axis.ticks.y = element_line(size = 1, color = "black"),  # Add y-axis ticks inside
        axis.ticks.length = unit(-0.3, "cm"),  # Move ticks inside the plot
        
        # Set a uniform font size for all text elements
        axis.title.x = element_text(size = 16),  # x-axis label
        axis.title.y = element_text(size = 16),  # y-axis label
        axis.text.x = element_text(size = 16),  # x-axis tick labels
        axis.text.y = element_text(size = 16),  # y-axis tick labels
        plot.title = element_text(size = 16, face = "bold"),  # Title font size and bold
        plot.margin = margin(1, 1, 1, 1, "cm")  # Add some margin around the plot for better spacing
      )
    
    # Subset Seurat object based on sample group and extract the feature data for "post" and "pre" samples
    post_values <- sapply(post_samples, function(current_sample) {
        # Subset based on the sample metadata and fetch the feature data
        subset_data <- subset(seurat_objects, subset = sample == current_sample)
        observation_of_interest <- as.matrix(FetchData(subset_data, vars = current_obs))
        spatial_coordinates <- SeuratObject::GetTissueCoordinates(subset_data)
        coordinates_matrix <- as.matrix(spatial_coordinates[, c('x', 'y')])
        moransI <- RunMoransI(t(observation_of_interest), coordinates_matrix, verbose = FALSE)
        moransI <- moransI$observed
        })
    
    pre_values <- sapply(pre_samples, function(current_sample) {
        subset_data <- subset(seurat_objects, subset = sample == current_sample)
        observation_of_interest <- as.matrix(FetchData(subset_data, vars = current_obs))
        spatial_coordinates <- SeuratObject::GetTissueCoordinates(subset_data)
        coordinates_matrix <- as.matrix(spatial_coordinates[, c('x', 'y')])
        moransI <- RunMoransI(t(observation_of_interest), coordinates_matrix, verbose = FALSE)
        moransI <- moransI$observed
    })
    
    # Perform a t-test
    t_test_result <- t.test(pre_values, post_values, var.equal = FALSE)
    cat("T-statistic:", round(t_test_result$statistic, 4), "p-value:", round(t_test_result$p.value, 4), "\n")
    
    # Prepare data for plotting
    data <- data.frame(
      Sample = rep(c('Pre', 'Post'), times = c(length(pre_values), length(post_values))),
      Value = c(pre_values, post_values)
    )
    
    # Create a plot with strip plot, horizontal mean line, no grid, ticks inside the plot, and uniform font size
    morans_plot <- ggplot(data, aes(x = Sample, y = Value)) +
      geom_jitter(width = 0.1, alpha = 0.7, color = 'skyblue') +  # Blue circles for individual samples
      stat_summary(fun = "mean", geom = "crossbar", color = "red", size = 1.5, width = 0.2, fatten = 1) + # Red horizontal line for means
      labs(
        title = paste('Pre vs Post', current_obs, 'MoransI\nT-test p-value:', round(t_test_result$p.value, 4)),
        x = 'Sample Group',
        y = paste('MoransI', current_obs)
      ) +
      theme_minimal() +  # Minimal theme (removes some default elements)
      theme(
        panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(size = 1, color = "black"),  # Add axis lines
        axis.ticks.y = element_line(size = 1, color = "black"),  # Add y-axis ticks inside
        axis.ticks.length = unit(-0.3, "cm"),  # Move ticks inside the plot
        
        # Set a uniform font size for all text elements
        axis.title.x = element_text(size = 16),  # x-axis label
        axis.title.y = element_text(size = 16),  # y-axis label
        axis.text.x = element_text(size = 16),  # x-axis tick labels
        axis.text.y = element_text(size = 16),  # y-axis tick labels
        plot.title = element_text(size = 16, face = "bold"),  # Title font size and bold
        plot.margin = margin(1, 1, 1, 1, "cm")  # Add some margin around the plot for better spacing
      )
    
    # Arrange the plots in one row
    grid_plot_lower <- grid.arrange(dist_plot, means_plot, morans_plot, ncol = 3)
    
    # Combine the upper and lower grid plots vertically with adjusted row heights
    grid_plot_combined <- grid.arrange(
      grid_plot_upper, 
      grid_plot_lower, 
      nrow = 2, 
      heights = c(3, 2)  # Adjust height of the rows (2 parts for upper, 1 part for lower)
    )
    
    output_path <- paste0(getwd(),'/figures/',current_metric,'_summary.png')
    
    ggsave(filename = output_path, plot = grid_plot_combined, width = 16, height = 12, dpi = 300)
    
}