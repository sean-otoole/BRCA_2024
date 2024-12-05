library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)

object_path <- paste0(getwd(),'/objects/','seurat_objects.rds')
seurat_objects <- readRDS(object_path)

# Create the individual VlnPlot objects
plot1 <- VlnPlot(seurat_objects, features = "nFeature_Spatial", group.by = 'sample', pt.size = 0)
plot2 <- VlnPlot(seurat_objects, features = "percent.mt", group.by = 'sample', pt.size = 0)
plot3 <- VlnPlot(seurat_objects, features = "nCount_Spatial", group.by = 'sample', pt.size = 0)

# Arrange the plots in a 1x3 grid layout
grid_plot <- grid.arrange(plot1, plot2, plot3, ncol = 3)

output_path <- paste0(getwd(),'/figures/','vlnplot_grid.png')

ggsave(output_path, plot = grid_plot, width = 16, height = 6, dpi = 300)