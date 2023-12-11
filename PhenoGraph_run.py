import phenograph
import numpy as np
import pandas as pd

data1_path = '/mnt/data1/fatema/IMC_T1D_data1/'
data2_path = '/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/'
data3_path= '/mnt/data1/fatema/IMC_T1Dsample/'
type = ['mgDF', 'mgDF_GAD+GAD-', 'T1D']


##############################################################################
type_id = 0

# data should be cell vs gene ndarray
data = np.load(data1_path + "gene_vs_cell_count_"+type[type_id]+".npy")
print(data.shape)
data = data.transpose()
print(data.shape)

df_normalized_protein = pd.read_csv('/mnt/data1/fatema/IMC_T1D_data1/normalized_protein.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins
data = df_normalized_protein.values
communities, graph, Q = phenograph.cluster(data) #higher k = less number of clusters. Default k=30


#####################################################################################

type_id = 2
# data should be cell vs gene ndarray
df_normalized_protein = pd.read_csv(data3_path+'normalized_protein_'+ type[type_id] +'.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins
data = df_normalized_protein.values
print(data.shape)
communities, graph, Q = phenograph.cluster(data) #higher k = less number of clusters. Default k=30

######################################################################################
type_id = 1
# data should be cell vs gene ndarray
df_normalized_protein = pd.read_csv(data2_path+'normalized_protein_'+ type[type_id] +'.csv', sep=",", header=0,index_col=0) 
df_normalized_protein = df_normalized_protein.T # now rows are cells and columns are proteins
data = df_normalized_protein.values
print(data.shape)
communities, graph, Q = phenograph.cluster(data, k=100) #higher k = less number of clusters. Default k=30
