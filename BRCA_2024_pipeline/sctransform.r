cat("R script started\n")

library(Seurat)
library(future)
options(future.globals.maxSize = 20 * 1024^3)  #15 GiB
plan("multicore", workers = 4)

objects_path <- paste0(getwd(), '/objects/seurat_objects.rds')

seurat_objects = readRDS(objects_path)
seurat_objects <- SCTransform(seurat_objects, assay = "Spatial", verbose = TRUE, vst.flavor = "v2", conserve.memory=FALSE,
                             return.only.var.genes = TRUE,  vars.to.regress = c('percent.mt'))  # originally return only variable genes was False
saveRDS(seurat_objects, objects_path)
cat("R script completed successfully\n")