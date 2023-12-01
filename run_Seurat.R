library(Seurat)             
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)

options = 'dt-randomCCC_equally_spaced_lrc105_cp100_noise0_threshold_dist_cellCount3000'

df=read.csv(file = paste("/cluster/home/t116508uhn/synthetic_cell_",options,"_x.csv",sep=""), header = FALSE) #read.csv(file = '/cluster/home/t116508uhn/synthetic_cell_type6_f_x.csv', header = FALSE)
cell_x=list()  
for(i in 1:ncol(df)) {      
  cell_x[[i]] <- df[ , i]    
}
df=read.csv(file = paste('/cluster/home/t116508uhn/synthetic_cell_',options,'_y.csv',sep=""), header = FALSE) #read.csv(file = '/cluster/home/t116508uhn/synthetic_cell_type6_f_y.csv', header = FALSE)
cell_y=list()  
for(i in 1:ncol(df)) {      
  cell_y[[i]] <- df[ , i]    
}

coord_x <- read.csv('/mnt/data1/fatema/x_coord_mgDF.csv', header = FALSE)
coord_x_vector <- coord_x[, 1]
coord_y <- read.csv('/mnt/data1/fatema/y_coord_mgDF.csv', header = FALSE)
coord_y_vector <- coord_y [, 1]






data1_path = '/mnt/data1/fatema/IMC_T1D_data1/'
countsData <- read.csv(file = paste('/mnt/data1/fatema/IMC_T1D_data1/gene_vs_cell_mgDF.csv', sep=""), row.names = 1) 

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

mat@meta.data$x <- coord_x_vector
mat@meta.data$y <- coord_y_vector


countsData <- mat
#####################################################################################
temp <- CreateSeuratObject(counts = countsData)
temp <- SCTransform(temp)
temp <- NormalizeData(temp)
temp <- FindVariableFeatures(temp) 
temp <- ScaleData(temp)

temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:30)
temp <- FindClusters(temp)
save(temp, file='/mnt/data1/fatema/IMC_T1D_data1/temp.Rda')
temp_matrix <- temp[['seurat_clusters']]
write.csv(temp_matrix, '/mnt/data1/fatema/IMC_T1D_data1/seurat_clusters.csv')

temp <- RunUMAP(temp , reduction = "pca", dims = 1:30)
p1 <- DimPlot(temp, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
p2 <- SpatialDimPlot(temp, label = TRUE,group.by = 'seurat_clusters', label.size = 3)
ggsave("/cluster/home/t116508uhn/64630/myplot.png", plot = (p1+p2))
