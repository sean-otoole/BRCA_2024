library(Seurat)
# library(future)
library(CellChat)
options(stringsAsFactors = FALSE)


setwd("/tmp/work/Visium/BRCA_2024")
options(future.globals.maxSize = 50 * 1024^3)  #50 GiB
# plan("multicore", workers = 30)

seurat_objects_list <- readRDS('seurat_objects_list.rds')

cellchat_list <- list()

dataset_paths = c('2018838/PRE-01','2018838/PRE-02','2018838/POST-03','2018838/POST-05',
           '2018839/PRE-04','2018839/PRE-05','2018839/POST-06','2018839/POST-11')

dataset_paths <- lapply(dataset_paths, function(x){paste('/tmp/work/Visium/', x, sep = '')})

for (seurat in seurat_objects_list) {

    #get the variable name ready for later
    var_name <- paste(unique(seurat[[]]$'sample'),'cellchat',sep='_')

    #construct the cellchat object
    
    #spatial coordinates
    spatial.locs = GetTissueCoordinates(seurat, scale = NULL, cols = c("imagerow", "imagecol"))
    #scale factors
    
    current_path <- dataset_paths[grep(paste0(unique(seurat$sample), "$"), dataset_paths)]
    current_path <- paste(current_path,'spatial',sep='/')

    scalefactors = jsonlite::fromJSON(txt = file.path(current_path, 'scalefactors_json.json'))
    spot.size = 65 # the theoretical spot size (um) in 10X Visium
    conversion.factor = spot.size/scalefactors$spot_diameter_fullres
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    d.spatial <- computeCellDistance(coordinates = as.matrix(spatial.locs[, c(1, 2)]), ratio = spatial.factors$ratio, tol = spatial.factors$tol)
    min(d.spatial[d.spatial!=0]) # this value should approximately equal 100um for 10X Visium data

    data.input = GetAssayData(seurat, slot = "data", assay = "SCT") # normalized data matrix
    meta = data.frame(labels = Idents(seurat), row.names = names(Idents(seurat))) # manually create a dataframe consisting of the cell labels
    
    if (any(unique(meta) == 0)) {
        meta$labels <- factor(as.numeric(meta$labels)+1)
    }

    #create the cellchat object
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = as.matrix(spatial.locs[,c(1,2)]), spatial.factors = spatial.factors)

    # set the database#
    cellchat@DB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    #communication probability in intracellular network
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.range = 200, scale.distance = 0.01, contact.range = 100, nboot = 20)

    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)

    #calculate the aggregated cell-cell communicaiton network
    cellchat <- aggregateNet(cellchat)
    
    assign(var_name, cellchat) 
    cellchat_list[[var_name]] <- get(var_name)

}

saveRDS(cellchat_list,'cellchat_list.rds')