##########
# Noah Burget
# 9/6/2023
# This script will determine the most DEP (diff. expressed protein) from a matrix w/ group/cluster labels
# Originally used for AnnoSpat paper analysis
##########
library(data.table)

proteinExpr.clusterLabeled <- read.table('/mnt/data1/fatema/IMC_T1D_data1/cluster_label_mgDF.csv',sep=',',header=T)[,-1] # File of cells x protein expression
clusters <- unique(proteinExpr.clusterLabeled$label)
test.results.foldChange <- data.frame()
test.results.wilcox <- data.frame()
for(i in clusters){ # separate cells in target cluster from all other cells
  target.cluster <- subset(proteinExpr.clusterLabeled, label == i)[,-c(1,35)] # Get cells in current target cluster
  other.clusters <- subset(proteinExpr.clusterLabeled, label != i)[,-c(1,35)] # Get cells in every other cluster
  for(j in colnames(target.cluster)){ # For every protein measured, compare 
    print(paste0("Testing ",j," in cluster ",i," against all other clusters"))
    target <- target.cluster[,j] # target cluster expression
    other <- other.clusters[,j] # all other cluster expression
    log2FC <- log2(mean(target) / mean(other)) # Calculate log2 fold change from mean expression values
    test.results.foldChange[i+1,j] <- log2FC
    res <- wilcox.test(target, other) # Use wilcox rank-sum test to determine p-value
    test.results.wilcox[i+1,j] <- res$p.value
  }
}

test.results.wilcox[test.results.wilcox == 0] <- .Machine$double.xmin # p-values = 0 are set to lowest double val
test.results.wilcox <- -log10(test.results.wilcox)
test.results.wilcoxFoldChange <- test.results.wilcox * test.results.foldChange
write.csv(test.results.wilcoxFoldChange,"/mnt/data1/fatema/IMC_T1D_data1/cluster_annotation_mgDF.csv.csv") #, row.names=FALSE
################################################################################################################################

proteinExpr.clusterLabeled <- read.table('/mnt/data1/fatema/IMC_T1Dsample/cluster_label_T1D.csv',sep=',',header=T)[,-1] # File of cells x protein expression
clusters <- unique(proteinExpr.clusterLabeled$label)
test.results.foldChange <- data.frame()
test.results.wilcox <- data.frame()
for(i in clusters){ # separate cells in target cluster from all other cells
  target.cluster <- subset(proteinExpr.clusterLabeled, label == i)[,-c(1,35)] # Get cells in current target cluster
  other.clusters <- subset(proteinExpr.clusterLabeled, label != i)[,-c(1,35)] # Get cells in every other cluster
  for(j in colnames(target.cluster)){ # For every protein measured, compare 
    print(paste0("Testing ",j," in cluster ",i," against all other clusters"))
    target <- target.cluster[,j] # target cluster expression
    other <- other.clusters[,j] # all other cluster expression
    log2FC <- log2(mean(target) / mean(other)) # Calculate log2 fold change from mean expression values
    test.results.foldChange[i+1,j] <- log2FC
    res <- wilcox.test(target, other) # Use wilcox rank-sum test to determine p-value
    test.results.wilcox[i+1,j] <- res$p.value
  }
}

test.results.wilcox[test.results.wilcox == 0] <- .Machine$double.xmin # p-values = 0 are set to lowest double val
test.results.wilcox <- -log10(test.results.wilcox)
test.results.wilcoxFoldChange <- test.results.wilcox * test.results.foldChange
#write.csv(test.results.wilcoxFoldChange,"/mnt/data1/fatema/IMC_T1Dsample/cluster_annotation_T1D.csv.csv") #, row.names=FALSE
temp <- transpose(test.results.wilcoxFoldChange)
colnames(temp) <- rownames(test.results.wilcoxFoldChange)
rownames(temp) <- colnames(test.results.wilcoxFoldChange)
write.csv(temp,"/mnt/data1/fatema/IMC_T1Dsample/cluster_annotation_T1D.csv") #, row.names=FALSE

####################################################################################################################################
proteinExpr.clusterLabeled <- read.table('/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/cluster_label_mgDF_GAD+GAD-.csv',sep=',',header=T)[,-1] # File of cells x protein expression
clusters <- unique(proteinExpr.clusterLabeled$label)
test.results.foldChange <- data.frame()
test.results.wilcox <- data.frame()
for(i in clusters){ # separate cells in target cluster from all other cells
  target.cluster <- subset(proteinExpr.clusterLabeled, label == i)[,-c(1,35)] # Get cells in current target cluster
  other.clusters <- subset(proteinExpr.clusterLabeled, label != i)[,-c(1,35)] # Get cells in every other cluster
  for(j in colnames(target.cluster)){ # For every protein measured, compare 
    print(paste0("Testing ",j," in cluster ",i," against all other clusters"))
    target <- target.cluster[,j] # target cluster expression
    other <- other.clusters[,j] # all other cluster expression
    log2FC <- log2(mean(target) / mean(other)) # Calculate log2 fold change from mean expression values
    test.results.foldChange[i+1,j] <- log2FC
    res <- wilcox.test(target, other) # Use wilcox rank-sum test to determine p-value
    test.results.wilcox[i+1,j] <- res$p.value
  }
}

test.results.wilcox[test.results.wilcox == 0] <- .Machine$double.xmin # p-values = 0 are set to lowest double val
test.results.wilcox <- -log10(test.results.wilcox)
test.results.wilcoxFoldChange <- test.results.wilcox * test.results.foldChange
write.csv(test.results.wilcoxFoldChange,"/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/cluster_annotation_mgDF_GAD+GAD-.csv") #, row.names=FALSE
