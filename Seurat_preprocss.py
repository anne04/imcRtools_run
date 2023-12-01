import os
import pandas as pd
import copy
import csv
import numpy as np
import sys
from collections import defaultdict
data1_path = '/mnt/data1/fatema/IMC_T1D_data1/'

df = pd.read_csv('/mnt/data1/fatema/IMC_T1D/raw_data/mgDF.csv', sep=",", header=0,index_col=0) 
# seperate the protein names, (x, y), and cell names

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
    cell_key = str(j+1)+'-'+file_name[j] + '-' + status_list[j]
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

gene_vs_cell = np.zeros((len(protein_name),len(cell_name)))
for i in range (0, len(protein_name)):
    gene = protein_name[i]
    for j in range (0, len(cell_name)):
        cell = cell_name[j]
        gene_vs_cell[i][j] = df[gene][cell]

###########
data_list=defaultdict(list)
gene_vs_cell = np.zeros((len(protein_name),len(cell_name)))
for i in range (0, len(protein_name)):
    gene = protein_name[i]
    for j in range (0, len(cell_name)):
        cell = cell_name[j]
        gene_vs_cell[i][j] = df[gene][cell]
        data_list[cell].append(df[gene][cell]) 
     
        
data_list_pd = pd.DataFrame(data_list)        
data_list_pd[' '] = protein_name
data_list_pd = data_list_pd.set_index(' ')    
data_list_pd.to_csv(data1_path+'gene_vs_cell_mgDF.csv')
np.save(data1_path+"gene_vs_cell_count_mgDF", gene_vs_cell)

df = pd.DataFrame(protein_name)
df.to_csv(data1_path+'protein_marker_mgDF.csv', index=False, header=False)

df = pd.DataFrame(cell_name)
df.to_csv(data1_path+'cell_id_mgDF.csv', index=False, header=False)

df = pd.DataFrame(file_name)
df.to_csv(data1_path+'file_name_mgDF.csv', index=False, header=False)

df = pd.DataFrame(status_list)
df.to_csv(data1_path+'status_list_mgDF.csv', index=False, header=False)

df = pd.DataFrame(x_coord)
df.to_csv(data1_path+'x_coord_mgDF.csv', index=False, header=False)

df = pd.DataFrame(y_coord)
df.to_csv(data1_path+'y_coord_mgDF.csv', index=False, header=False)

df = pd.DataFrame(cell_label)
df.to_csv(data1_path+'cell_label_islets_mgDF.csv', index=False, header=False)

