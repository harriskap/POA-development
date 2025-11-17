library(Seurat)
library(Augur)
library(glue)

name_data <- "inhibitory_merged_noNA.RDS"
name_out <- "default_param"

# Subsetting dataset by age
args <- commandArgs(trailingOnly=TRUE)
age_ <- args[1]
print(age_)

# Load data
seurat <- readRDS(glue("../{name_data}"))
DefaultAssay(seurat) <- "RNA"

# Remove sex related genes
rm_genes <- read.csv("../y_chromosome_and_x_inactivation_genes.csv", header=FALSE)
rm_genes <- rm_genes$V1
seurat <- seurat[setdiff(x = rownames(x = seurat), y = rm_genes), ]

genes_mt <- rownames(seurat)[grepl("^mt-", rownames(seurat))]
genes_rb <- rownames(seurat)[grepl("^Rps|Rrpl", rownames(seurat))]

seurat <- seurat[setdiff(x = rownames(x = seurat), y = genes_mt), ]
seurat <- seurat[setdiff(x = rownames(x = seurat), y = genes_rb), ]

# Subset seurat object to one age
# Re-normalize
seurat_sub <- subset(seurat, subset = (age == age_))
seurat_sub <- NormalizeData(seurat_sub, assay="RNA")
seurat_sub <- ScaleData(object = seurat_sub, assay="RNA", features=rownames(seurat_sub))

# Run augur
#augur = calculate_auc(seurat_sub, seurat_sub@meta.data, cell_type_col = "my.cell.type", label_col = "sex", var_quantile=0.9, subsample_size=10,
#    n_threads = 16, rf_params = list(trees = 500, mtry = 2, min_n = NULL, importance = "accuracy"))
augur = calculate_auc(seurat_sub, seurat_sub@meta.data, cell_type_col = "my.cell.type", label_col = "sex", var_quantile=0.9, subsample_size=10,
    n_threads = 16)


saveRDS(augur, glue("{name_out}/{age_}.RDS"))

# Save AUC scores
df <- data.frame(augur$AUC)
colnames(df) <- c("cell_type", age_)
write.csv(glue("{name_out}/{age_}.csv"))
