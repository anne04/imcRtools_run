library(Seurat)             
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(reticulate)
np <- import("numpy")

####################################################################################
mat <- np$load("/mnt/data1/fatema/IMC_T1D_data1/gene_vs_cell_count_mgDF.npy")
cell_barcodes <- read.csv('/mnt/data1/fatema/IMC_T1D_data1/cell_id_mgDF.csv', header = FALSE)
# nrow(cell_barcodes) is 1406 and ncol(cell_barcodes) is 1.
cell_vector <- cell_barcodes[, 1]
gene_ids <- read.csv('/mnt/data1/fatema/IMC_T1D_data1/protein_marker_mgDF.csv', header = FALSE)
# nrow(gene_ids) is 19523 and ncol(gene_ids) is 1.
gene_vector <- gene_ids[, 1]
colnames(mat) <- cell_vector
rownames(mat) <- gene_vector
#mat@meta.data$x <- coord_x_vector
#mat@meta.data$y <- coord_y_vector
countsData <- mat
#####################################################################################
temp <- CreateSeuratObject(counts = countsData)
temp <- SCTransform(temp)
temp <- NormalizeData(temp)
temp <- FindVariableFeatures(temp) 
temp <- ScaleData(temp)

temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:30)
temp <- FindClusters(temp, resolution=0.3)
save(temp, file='/mnt/data1/fatema/IMC_T1D_data1/temp.Rda')
temp_matrix <- temp[['seurat_clusters']]
write.csv(temp_matrix, '/mnt/data1/fatema/IMC_T1D_data1/seurat_clusters.csv')
temp_protein_cell <- GetAssayData(object = temp, layer = "scale.data")
temp_protein_cell <- as.matrix(temp_protein_cell)
write.csv(temp_protein_cell, paste('/mnt/data1/fatema/IMC_T1D_data1/scaled_protein.csv',sep=""))

temp_protein_cell <- GetAssayData(object = temp, layer = "data")
temp_protein_cell <- as.matrix(temp_protein_cell)
write.csv(temp_protein_cell, paste('/mnt/data1/fatema/IMC_T1D_data1/normalized_protein.csv',sep=""))

#################################################################################

####################################################################################
mat <- np$load("/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/gene_vs_cell_count_mgDF_GAD+GAD-.npy")
cell_barcodes <- read.csv('/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/cell_id_mgDF_GAD+GAD-.csv', header = FALSE)
# nrow(cell_barcodes) is 1406 and ncol(cell_barcodes) is 1.
cell_vector <- cell_barcodes[, 1]
gene_ids <- read.csv('/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/protein_marker_mgDF_GAD+GAD-.csv', header = FALSE)
# nrow(gene_ids) is 19523 and ncol(gene_ids) is 1.
gene_vector <- gene_ids[, 1]
colnames(mat) <- cell_vector
rownames(mat) <- gene_vector
#mat@meta.data$x <- coord_x_vector
#mat@meta.data$y <- coord_y_vector
countsData <- mat
#####################################################################################
temp <- CreateSeuratObject(counts = countsData)
#temp <- SCTransform(temp)
temp <- NormalizeData(temp) #https://htmlpreview.github.io/?https://github.com/satijalab/sctransform/blob/supp_html/supplement/seurat.html
temp <- FindVariableFeatures(temp) 
temp <- ScaleData(temp)

temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:30)
temp <- FindClusters(temp, resolution=0.3)
save(temp, file='/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/mgDF_GAD+GAD-temp.Rda')
temp_matrix <- temp[['seurat_clusters']]
write.csv(temp_matrix, '/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/seurat_clusters_mgDF_GAD+GAD-.csv')
#################################################################################


####################################################################################
mat <- np$load("/mnt/data1/fatema/IMC_T1Dsample/gene_vs_cell_count_T1D.npy")
cell_barcodes <- read.csv('/mnt/data1/fatema/IMC_T1Dsample/cell_id_T1D.csv', header = FALSE)
# nrow(cell_barcodes) is 1406 and ncol(cell_barcodes) is 1.
cell_vector <- cell_barcodes[, 1]
gene_ids <- read.csv('/mnt/data1/fatema/IMC_T1Dsample/protein_marker_T1D.csv', header = FALSE)
# nrow(gene_ids) is 19523 and ncol(gene_ids) is 1.
gene_vector <- gene_ids[, 1]
colnames(mat) <- cell_vector
rownames(mat) <- gene_vector
#mat@meta.data$x <- coord_x_vector
#mat@meta.data$y <- coord_y_vector
countsData <- mat
#####################################################################################
temp <- CreateSeuratObject(counts = countsData)
#temp <- SCTransform(temp)
temp <- NormalizeData(temp) #https://htmlpreview.github.io/?https://github.com/satijalab/sctransform/blob/supp_html/supplement/seurat.html
temp <- FindVariableFeatures(temp) 
temp <- ScaleData(temp)

temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:30)
temp <- FindClusters(temp, resolution=0.3)
save(temp, file='/mnt/data1/fatema/IMC_T1Dsample/T1D-temp.Rda')
temp_matrix <- temp[['seurat_clusters']]
write.csv(temp_matrix, '/mnt/data1/fatema/IMC_T1Dsample/seurat_clusters_T1D.csv')
#################################################################################










temp <- RunUMAP(temp , reduction = "pca", dims = 1:30)
p1 <- DimPlot(temp, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
p2 <- SpatialDimPlot(temp, label = TRUE,group.by = 'seurat_clusters', label.size = 3)
ggsave("/cluster/home/t116508uhn/64630/myplot.png", plot = (p1+p2))
