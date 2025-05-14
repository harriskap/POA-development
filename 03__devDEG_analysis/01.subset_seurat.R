# export data for eigentrends
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
# https://github.com/mojaveazure/seurat-disk/issues/102

library(SeuratDisk)
library(Seurat)
#library(SeuratObject)
library(tidyverse)

# prepare excitatory h5 object

# excitatory, inhibitory, merged

experiment <- "merged"

if (experiment == "excitatory"){
    seurat_input_file <- "data/2023-05-16_input/excitatory_merged_noNA.RDS"
    seurat_subset_file <- "data/2023-05-16_input/seurat.excitatory_select.RDS"
    seurat_h5_file <- "data/2023-05-16_input/seurat.excitatory_select.h5Seurat"
    subset_cell_types <- c("e-A1", "e-A2", "e-A3",
                           "e-F1",
                           "e-M1", "e-M2", "e-M3", "e-M4", 
                           "e-M5", "e-M8", "e-M9", "e-M10",
                           "e-P1", "e-P2", "e-P3")
}else if (experiment == "inhibitory"){
    seurat <- "data/2023-05-16_input/inhibitory_merged_noNA.RDS"
    seurat_subset_file <- "data/2023-05-16_input/seurat.inhibitory_select.RDS"
    seurat_h5_file <- "data/2023-05-16_input/seurat.inhibitory_select.h5Seurat"
    # removing i-M8, i-S2 due to low pseudobulk numbers
    subset_cell_types <- c("i-B1", "i-B2", "i-B3",
                           "i-H1", "i-H2", "i-H3",
                           "i-M1", "i-M2", "i-M3",  
                           "i-M4", "i-M5", "i-M6",
                           "i-M7", "i-M9", 
                           "i-P1", 
                           "i-S1", "i-S3", "i-S4",
                           "i-S5",
                           "i-V1")
}else{
   seurat_subset_file <- "data/2023-05-16_input/seurat.merged.RDS"
   seurat_h5_file <- "data/seurat.merged.h5Seurat"
}


if (experiment %in% c("excitatory", "inhibitory")){
    if (file.exists(seurat_subset_file)){
        seurat_subset <- readRDS(seurat_subset_file)
    }else{
        seurat <- readRDS(seurat_input_file)
        # first we selected by region
        # seurat_avpe <- subset(seurat, subset = brain.region == "AvPE.MnPO")
        # selecting by cell types 
        seurat_subset <- subset(seurat, subset = my.cell.type %in% subset_cell_types)
        saveRDS(seurat_subset, seurat_subset_file)
    }
}else{
    # this needs more RAM than laptop
    # 50G was not enough to merge
    # 200G
    if (file.exists(seurat_subset_file)){
        seurat_subset <- readRDS(seurat_subset_file)
    }else{
        seurat <- readRDS("data/2023-05-16_input/excitatory_merged_noNA.RDS")
        seurat2 <- readRDS("data/2023-05-16_input/inhibitory_merged_noNA.RDS")
        seurat <- merge(seurat, y = seurat2, merge.data = TRUE)
    }
    seurat_subset <- seurat
}
  
if (!file.exists(seurat_h5_file)){
      slot(seurat$SCT@SCTModel.list[[1]], 'median_umi') = median(seurat$SCT@SCTModel.list[[1]]@cell.attributes$umi)
      DefaultAssay(seurat_subset) <- "RNA"
      SaveH5Seurat(seurat_subset, seurat_h5_file)
      Convert(seurat_h5_file,
          dest = "h5ad", 
          assay = "RNA")
}
# write10xCounts(x = seurat_subset@assays$RNA@counts, "data/dev_deg_v4/excitatory_counts")


# was needed for hotspot method
# write10xCounts(x = seurat_subset@assays$RNA@counts, "data/dev_deg_v4/excitatory_counts"


# DimPlot(seurat_subset)
# mdata <- seurat_avpe@meta.data