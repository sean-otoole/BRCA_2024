library(Seurat)
library(dplyr)
library(future)
setwd("/tmp/work/Visium/BRCA_2024/objects")
options(future.globals.maxSize = 20 * 1024^3)  #50 GiB
plan("multicore", workers = 4)
ref_data = readRDS('/tmp/work/Visium/zenodo.11468564/scRNA_data/panB_scRNA_processed_data.rds')
ref_data <- SCTransform(ref_data, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
saveRDS('b_cell_ref_data.rds')