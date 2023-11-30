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

countsData <- read.csv(file = paste('/cluster/home/t116508uhn/synthetic_gene_vs_cell_',options,'.csv', sep=""),row.names = 1) # read.csv(file = '/cluster/home/t116508uhn/synthetic_gene_vs_cell_type6_f.csv',row.names = 1)
pdac_sample <- CreateSeuratObject(counts = countsData)
#temp <- SCTransform(pdac_sample)
temp <- ScaleData(pdac_sample)
temp <- FindVariableFeatures(temp) 
temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:30)
temp <- FindClusters(temp, verbose = FALSE)
temp <- RunUMAP(temp , reduction = "pca", dims = 1:30)
p1 <- DimPlot(temp, reduction = "umap",group.by = 'seurat_clusters', label = TRUE)
p2 <- SpatialDimPlot(temp, label = TRUE,group.by = 'seurat_clusters', label.size = 3)
ggsave("/cluster/home/t116508uhn/64630/myplot.png", plot = (p1+p2))
