import os
import pandas as pd
import copy
import csv
import numpy as np
import sys
from collections import defaultdict

data1_path_to = '/mnt/data1/fatema/IMC_T1D_data1/'
data2_path_to = '/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/'
data3_path_to = '/mnt/data1/fatema/IMC_T1Dsample/'

type = ['mgDF', 'mgDF_GAD+GAD-', 'T1D']


type_id = 0
df_normalized_protein = pd.read_csv(data1_path_to+'normalized_protein.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins

clusters = np.load(data1_path_to+'phenograph_clusters_less_'+type[type_id]+'.npy')  # ? clusters with k = 1000
num_cells = clusters.shape[0]
cluster_id = []
for i in range (0, num_cells):
    cluster_id.append(clusters[i])


df = pd.read_csv('/mnt/data1/fatema/seurat_clusters/mgDF_data/seurat_clusters_mgDF.csv', sep=",", header=0,index_col=0) # 22 clusters with resolution = 0.2
df['seurat_clusters']=cluster_id
df = df.rename(columns={'seurat_clusters': 'phenograph_clusters'})
df.to_csv('/mnt/data1/fatema/phenograph_clusters/mgDF_data/'+'phenograph_clusters_'+type[type_id]+'.csv')


df_normalized_protein['label']=cluster_id   
df_normalized_protein.to_csv(data1_path_to+'phenograph_cluster_label_'+type[type_id]+'.csv')

############################################
type_id = 1
df_normalized_protein = pd.read_csv(data2_path_to+'normalized_protein.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins

clusters = np.load(data2_path_to+'phenograph_clusters_less_'+type[type_id]+'.npy') # 23 with k =1000.
num_cells = clusters.shape[0]
cluster_id = []
for i in range (0, num_cells):
    cluster_id.append(clusters[i])

df = pd.read_csv('/mnt/data1/fatema/seurat_clusters/mgDF_GAD+-_data/seurat_clusters_mgDF_GAD+GAD-.csv', sep=",", header=0,index_col=0) # 19 clusters with resolution-0.2
df['seurat_clusters']=cluster_id
df = df.rename(columns={'seurat_clusters': 'phenograph_clusters'})
df.to_csv('/mnt/data1/fatema/phenograph_clusters/mgDF_GAD+-_data/'+'phenograph_clusters_'+type[type_id]+'.csv')

df_normalized_protein['label']=cluster_id   
df_normalized_protein.to_csv(data2_path_to+'phenograph_cluster_label_'+type[type_id]+'.csv')
######################################################################################

type_id = 2
df_normalized_protein = pd.read_csv(data3_path_to+'normalized_protein.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins

clusters = np.load(data3_path_to+'phenograph_clusters_less_'+type[type_id]+'.npy')  # 21 clusters with k = 1000
num_cells = clusters.shape[0]
cluster_id = []
for i in range (0, num_cells):
    cluster_id.append(clusters[i])


df = pd.read_csv('/mnt/data1/fatema/seurat_clusters/T1D_data/seurat_clusters_T1D.csv', sep=",", header=0,index_col=0) # 17 clusters with resolution = 0.3
df['seurat_clusters']=cluster_id
df = df.rename(columns={'seurat_clusters': 'phenograph_clusters'})
df.to_csv('/mnt/data1/fatema/phenograph_clusters/T1D_data/'+'phenograph_clusters_'+type[type_id]+'.csv')


df_normalized_protein['label']=cluster_id   
df_normalized_protein.to_csv(data3_path_to+'phenograph_cluster_label_'+type[type_id]+'.csv')
