import phenograph
import numpy as np


data1_path = '/mnt/data1/fatema/IMC_T1D_data1/'
data2_path = '/mnt/data1/fatema/IMC_mgdf_gad+_gad-sample/'
type = ['mgDF', 'mgDF_GAD+GAD-', 'T1D']
type_id = 0

#data1_path_from = '/mnt/data1/fatema/IMC_T1D/raw_data/mgDF.csv'
#data2_path_from = '/mnt/data1/fatema/IMC_T1D/raw_data/mgDF_GAD+GAD-samples.csv'



# data should be cell vs gene ndarray
data = np.load(data1_path + "gene_vs_cell_count_"+type[type_id], gene_vs_cell)
print(data.shape)
data.transpose()
print(data.shape)

communities, graph, Q = phenograph.cluster(data)
