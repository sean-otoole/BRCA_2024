library(Seurat)
library(dplyr)
library(future)
setwd("/tmp/work/Visium/BRCA_2024/objects")
options(future.globals.maxSize = 130 * 1024^3)  #15 GiB
plan("multicore", workers = 4)
ref_data = readRDS('/tmp/work/Visium/zenodo.11468564/scRNA_data/panB_scRNA_processed_data.rds')
ref_data <- SCTransform(ref_data, ncells = 3000, verbose = TRUE, conserve.memory=TRUE) %>%
    RunPCA(verbose = TRUE) %>%
    RunUMAP(dims = 1:30)

saveRDS(ref_data,'b_cell_ref_data.rds')

visium <- readRDS('/tmp/work/Visium/BRCA_2024/objects/seurat_objects.rds')
visium <- RunPCA(object = visium)

visium_split <- SplitObject(visium, split.by = "sample")

for (i in seq_along(visium_split)) {
    item <- visium_split[[i]]
    anchors <- FindTransferAnchors(reference = ref_data, query = item, normalization.method = "SCT", recompute.residuals = FALSE)
    predictions.assay <- TransferData(anchorset = anchors, refdata = ref_data$celltype, 
                                      prediction.assay = TRUE, weight.reduction = item[["pca"]], dims = 1:30)
    item[["predictions"]] <- predictions.assay
    visium_split[[i]] <- item  # Store the modified object back into the list
}

visium <- merge(x = visium_split[[1]], y = visium_split[-1])

saveRDS(visium,'visium_labeled.rds')