library(reticulate)
library(SingleCellExperiment)
np <- import("numpy")
library(imcRtools)


mat <- np$load("/mnt/data1/fatema/gene_vs_cell_count_mgDF.npy")

cell_barcodes <- read.csv('/mnt/data1/fatema/cell_id_mgDF.csv', header = FALSE)
# nrow(cell_barcodes) is 1406 and ncol(cell_barcodes) is 1.
cell_vector <- cell_barcodes[, 1]


gene_ids <- read.csv('/mnt/data1/fatema/protein_marker_mgDF.csv', header = FALSE)
# nrow(gene_ids) is 19523 and ncol(gene_ids) is 1.
gene_vector <- gene_ids[, 1]


image_ids <- read.csv('/mnt/data1/fatema/file_name_mgDF.csv', header = FALSE)
image_id_vector <- image_ids[, 1]


coord_x <- read.csv('/mnt/data1/fatema/x_coord_mgDF.csv', header = FALSE)
coord_x_vector <- coord_x[, 1]

coord_y <- read.csv('/mnt/data1/fatema/y_coord_mgDF.csv', header = FALSE)
coord_y_vector <- coord_y [, 1]

cell_label <- read.csv('/mnt/data1/fatema/cell_label_mgDF.csv', header = FALSE)
cell_label_vector <- cell_label[, 1]

colnames(mat) <- cell_vector
rownames(mat) <- gene_vector

sce <- SingleCellExperiment(mat, colData=DataFrame(ImageNb=image_id_vector, Pos_X = coord_x_vector, Pos_Y = coord_y_vector, CellType = cell_label_vector))

pancreasSCE <- buildSpatialGraph(sce, img_id = "ImageNb", type = "knn", k = 3)

out <- testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", method = "histocat", colPairName = "knn_interaction_graph", iter = 1000)) #, BPPARAM = SerialParam(RNGseed = 123

save(out, file='/mnt/data1/fatema/out.Rda')

load(file='/mnt/data1/fatema/out.Rda')
write.csv(out, "/mnt/data1/fatema/out.csv", row.names=FALSE

############################################################
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
pretend.cell.labels <- sample(letters, ncol(counts), replace=TRUE)
pretend.gene.lengths <- sample(10000, nrow(counts))

sce2 <- SingleCellExperiment(list(counts=counts),
    colData=DataFrame(label=pretend.cell.labels),
    rowData=DataFrame(length=pretend.gene.lengths, length2=pretend.gene.lengths),
    metadata=list(study="GSE111111")
)


