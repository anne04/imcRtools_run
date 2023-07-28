import os
#import glob
import pandas as pd
#import shutil
import copy
import csv
import numpy as np
import sys

df_aanchal_label = pd.read_csv('/mnt/data1/gw/research/sc_integration/labels/labels_aanchal_imc_updated_1.csv', sep=",", header=0,index_col=0) 
label_dictionaly = dict()
for i in range (0, len(df_aanchal_label.index)):
    label_dictionaly[df_aanchal_label.index[i]] = df_aanchal_label['label'][df_aanchal_label.index[i]]


df = pd.read_csv('/mnt/data1/gw/research/sc_integration/data/imc/raw_data/mgDF.csv', sep=",", header=0,index_col=0) 
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
    if cell_key in label_dictionaly:
        cell_label.append(label_dictionaly[cell_key])
    else:
        cell_label.append('not_defined')
        not_defined_cell_count = not_defined_cell_count + 1

gene_vs_cell = np.zeros((len(protein_name),len(cell_name)))
for i in range (0, len(protein_name)):
    gene = protein_name[i]
    for j in range (0, len(cell_name)):
        cell = cell_name[j]
        gene_vs_cell[i][j] = df[gene][cell]

np.save("mnt/data1/fatema/gene_vs_cell_count_mgDF", gene_vs_cell)

df = pd.DataFrame(protein_name)
df.to_csv('mnt/data1/fatema/protein_marker_mgDF.csv', index=False, header=False)

df = pd.DataFrame(cell_name)
df.to_csv('mnt/data1/fatema/cell_id_mgDF.csv', index=False, header=False)

df = pd.DataFrame(file_name)
df.to_csv('mnt/data1/fatema/file_name_mgDF.csv', index=False, header=False)

df = pd.DataFrame(status_list)
df.to_csv('mnt/data1/fatema/status_list_mgDF.csv', index=False, header=False)

df = pd.DataFrame(x_coord)
df.to_csv('mnt/data1/fatema/x_coord_mgDF.csv', index=False, header=False)

df = pd.DataFrame(y_coord)
df.to_csv('mnt/data1/fatema/y_coord_mgDF.csv', index=False, header=False)


df = pd.DataFrame(cell_label)
df.to_csv('/mnt/data1/fatema/cell_label_mgDF.csv', index=False, header=False)
   
