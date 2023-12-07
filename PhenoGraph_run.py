import phenograph
import numpy as np
import pandas as pd

data1_path = '/mnt/data1/fatema/IMC_T1D_data1/'
data2_path = '/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/'
type = ['mgDF', 'mgDF_GAD+GAD-', 'T1D']
type_id = 0

# data should be cell vs gene ndarray
data = np.load(data1_path + "gene_vs_cell_count_"+type[type_id]+".npy")
print(data.shape)
data = data.transpose()
print(data.shape)

df_normalized_protein = pd.read_csv('/mnt/data1/fatema/IMC_T1D_data1/normalized_protein.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins

# remove first row and first column
df_normalized_protein.drop(index=df.index[0], axis=0, inplace=True)
df_normalized_protein.pop('column-name')
data = df_normalized_protein.to_numpy()


communities, graph, Q = phenograph.cluster(data) #higher k = less number of clusters. Default k=30
