library(Seurat)
library(future)
setwd("/tmp/work/Visium/BRCA_2024")
options(future.globals.maxSize = 60 * 1024^3)  # 10 GB
plan("multicore", workers = 50)

seurat_objects = readRDS('seurat_objects.rds')

seurat_objects_list <- SplitObject(seurat_objects, split.by = "sample")  #split the samples

seurat_objects_list <- lapply(seurat_objects_list, function(x) {
  if ("SCT" %in% Assays(x)) {
    FindSpatiallyVariableFeatures(x, assay = "SCT",  features = VariableFeatures(x)[1:3000],selection.method = "moransi")
  } else {
    warning("SCT assay not found in object")
    return(x)
  }
})

saveRDS(seurat_objects_list,'/tmp/work/Visium/BRCA_2024/seurat_objects_list.rds')