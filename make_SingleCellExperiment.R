library(reticulate)
library(SingleCellExperiment)
np <- import("numpy")

mat <- np$load("/mnt/data1/fatema/gene_vs_cell_count_mgDF.npy")

cell_barcodes <- read.csv('/mnt/data1/fatema/cell_id_mgDF.csv', header = FALSE)
# nrow(cell_barcodes) is 1406 and ncol(cell_barcodes) is 1.

gene_ids <- read.csv('/mnt/data1/fatema/protein_marker_mgDF.csv', header = FALSE)
# nrow(gene_ids) is 19523 and ncol(gene_ids) is 1.


cell_vector <- cell_barcodes[, 1]
gene_vector <- gene_ids[, 1]

colnames(mat) <- cell_vector
rownames(mat) <- gene_vector

sce <- SingleCellExperiment(counts)
