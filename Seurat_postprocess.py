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

df_clusters = pd.read_csv(data1_path_to+'seurat_clusters.csv', sep=",", header=0,index_col=0) 
num_cells = len(df_clusters.index)
cluster_id = []
for i in range (0, num_cells):
    cluster_id.append(df_clusters['seurat_clusters'][df_clusters.index[i]])
    
df_normalized_protein['label']=cluster_id   
df_normalized_protein.to_csv(data1_path_to+'cluster_label_'+type[type_id]+'.csv')

############################################
type_id = 1
df_normalized_protein = pd.read_csv(data2_path_to+'normalized_protein_'+type[type_id]+'.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins

df_clusters = pd.read_csv(data2_path_to+'seurat_clusters_'+ type[type_id] +'.csv', sep=",", header=0,index_col=0) 
num_cells = len(df_clusters.index)
cluster_id = []
for i in range (0, num_cells):
    cluster_id.append(df_clusters['seurat_clusters'][df_clusters.index[i]])
    
df_normalized_protein['label']=cluster_id   
df_normalized_protein.to_csv(data2_path_to+'cluster_label_'+type[type_id]+'.csv')

######################################################################################

type_id = 2
df_normalized_protein = pd.read_csv(data3_path_to+'normalized_protein_'+type[type_id]+'.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins

df_clusters = pd.read_csv(data3_path_to+'seurat_clusters_'+ type[type_id] +'.csv', sep=",", header=0,index_col=0) 
num_cells = len(df_clusters.index)
cluster_id = []
for i in range (0, num_cells):
    cluster_id.append(df_clusters['seurat_clusters'][df_clusters.index[i]])
    
df_normalized_protein['label']=cluster_id   
df_normalized_protein.to_csv(data3_path_to+'cluster_label_'+type[type_id]+'.csv')
