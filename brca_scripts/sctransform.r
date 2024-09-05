library(Seurat)
library(future)
setwd("/tmp/work/Visium/BRCA_2024")
options(future.globals.maxSize = 15 * 1024^3)  #15 GiB
plan("multicore", workers = 20)
seurat_objects = readRDS('seurat_objects.rds')
seurat_objects <- SCTransform(seurat_objects, assay = "Spatial", verbose = TRUE, vst.flavor = "v2", conserve.memory=FALSE)
saveRDS(seurat_objects,'seurat_objects.rds')