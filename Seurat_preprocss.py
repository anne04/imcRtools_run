import os
import pandas as pd
import copy
import csv
import numpy as np
import sys
from collections import defaultdict

data1_path_to = '/mnt/data1/fatema/IMC_T1D_data1/'
data2_path_to = '/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/'

data1_path_from = '/mnt/data1/fatema/IMC_T1D/raw_data/mgDF.csv'
data2_path_from = '/mnt/data1/fatema/IMC_T1D/raw_data/mgDF_GAD+GAD-samples.csv'

type = ['mgDF', 'mgDF_GAD+GAD-', 'T1D']
type_id = 0

df = pd.read_csv(data2_path_from, sep=",", header=0,index_col=0) 
# seperate the protein names, (x, y), and cell names

###################################
#df_column_names = list(df.columns)
#df.pop('column-name')

df_normalized_protein = pd.read_csv('/mnt/data1/fatema/IMC_T1D_data1/normalized_protein.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins

df_clusters = pd.read_csv('/mnt/data1/fatema/IMC_T1D_data1/seurat_clusters.csv', sep=",", header=0,index_col=0) 
num_cells = len(df_clusters.index)
cluster_id = []
for i in range (0, num_cells):
    cluster_id.append(df_clusters['seurat_clusters'][df_clusters.index[i]])
    
df_normalized_protein['label']=cluster_id   
df_normalized_protein.to_csv(data2_path_to+'cluster_label_'+type[type_id]+'.csv', index=False, header=False)

###########################
cell_name = []
x_coord = []
y_coord = []
for i in range (0, len(df.index)):
    cell_name.append(df.index[i])
    x_coord.append(df['CenterX'][df.index[i]])
    y_coord.append(df['CenterY'][df.index[i]])
    

protein_name = []
for i in range (2, len(df.columns)-2):
    protein_name.append(df.columns[i])

file_name = []
status_list = []
cell_label = []
not_defined_cell_count = 0
for j in range (0, len(cell_name)):
    cell = cell_name[j]
    file_name.append(df['TIFFfilename'][cell])
    status_list.append(df['Status'][cell])
'''    cell_key = str(j+1)+'-'+file_name[j] + '-' + status_list[j]
    if cell_key in label_dictionary:   
        # mark all the greek letter cell types as 'Islet_Cells'
        label = ''
        if label_dictionary[cell_key] in islet:
            label = 'Islet_Cells'
        else:
            label = label_dictionary[cell_key]
                   
        cell_label.append(label)
    else:
        cell_label.append('not_defined')
        not_defined_cell_count = not_defined_cell_count + 1

print('unique cell label: %d'%len(set(cell_label)))
print(set(cell_label))
'''

###########
data_list=defaultdict(list)
gene_vs_cell = np.zeros((len(protein_name),len(cell_name)))
for i in range (0, len(protein_name)):
    gene = protein_name[i]
    print(gene)
    for j in range (0, len(cell_name)):
        cell = cell_name[j]
        gene_vs_cell[i][j] = df[gene][cell]
        data_list[cell].append(df[gene][cell]) 
     
        
data_list_pd = pd.DataFrame(data_list)        
data_list_pd[' '] = protein_name
data_list_pd = data_list_pd.set_index(' ')    
data_list_pd.to_csv(data2_path_to+'gene_vs_cell_'+type[type_id]+'.csv')
np.save(data2_path_to+"gene_vs_cell_count_'+type[type_id]+'", gene_vs_cell)

df = pd.DataFrame(protein_name)
df.to_csv(data2_path_to+'protein_marker_'+type[type_id]+'.csv', index=False, header=False)

df = pd.DataFrame(cell_name)
df.to_csv(data2_path_to+'cell_id_'+type[type_id]+'.csv', index=False, header=False)

df = pd.DataFrame(file_name)
df.to_csv(data2_path_to+'file_name_'+type[type_id]+'.csv', index=False, header=False)

df = pd.DataFrame(status_list)
df.to_csv(data2_path_to+'status_list_'+type[type_id]+'.csv', index=False, header=False)

df = pd.DataFrame(x_coord)
df.to_csv(data2_path_to+'x_coord_'+type[type_id]+'.csv', index=False, header=False)

df = pd.DataFrame(y_coord)
df.to_csv(data2_path_to+'y_coord_'+type[type_id]+'.csv', index=False, header=False)

print('all done')
#df = pd.DataFrame(cell_label)
#df.to_csv(data2_path_to+'cell_label_islets_'+type[type_id]+'.csv', index=False, header=False)

