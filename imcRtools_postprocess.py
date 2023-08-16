import os
import pandas as pd
import copy
import csv
import numpy as np
import sys

save_path = '/mnt/data1/fatema/'
###################### read cell labels, status, ROI filenames #######################################################

cell_label=[]
with open('/mnt/data1/fatema/cell_label_islets_mgDF.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        cell_label.append(line[0])
        
status_list=[]
with open('/mnt/data1/fatema/status_list_mgDF.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        status_list.append(line[0])

file_name=[]
with open('/mnt/data1/fatema/file_name_mgDF.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        file_name.append(line[0])

######### find the list of regions for 'control' or 'GAD-' status, 'GAD+' status, '' #######################
ROI_control = dict()
ROI_AAB_poz = dict()
for cell_serial in range (0, len(file_name)):
    if status_list[cell_serial] == 'GAD-':
        ROI_control[file_name[cell_serial]] = ''
    elif status_list[cell_serial] == 'GAD+':
        ROI_AAB_poz[file_name[cell_serial]] = ''

print(ROI_control.keys())
print('count of ROI with control state is %d'%len(ROI_control.keys()))
ROI_control = list(ROI_control.keys())

########################### read histocat analysis result #############################################

out_histocat=[]
with open('/mnt/data1/fatema/out.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        out_histocat.append(line)

############### Islet vs Islet: find the 'ct' (and 'P'?) values for each of these ROI_control ########
from collections import defaultdict

ROI_control_Islet_Islet_ct = dict() # column 3
ROI_control_ct_distribution = defaultdict(list)
file_name_list = []
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Islet_Cells':
            ROI_control_Islet_Islet_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] == '1':
                ROI_control_ct_distribution['Distribution_type'].append('Islet vs Islet cells in Control')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
            elif out_histocat[i][9] == 'NA':
                continue
            else:
                ROI_control_ct_distribution['Distribution_type'].append('Islet vs Islet cells in Control')
                ROI_control_ct_distribution['ct'].append(0)
            file_name_list.append(out_histocat[i][0])
#################### draw density plot of ct values for islet vs islet cells in control ###############

df = pd.DataFrame(ROI_control_ct_distribution)
chart =alt.Chart(df).transform_density(
    'ct',
    as_=['ct', 'density'],
    groupby=['Distribution_type'],
    #bandwidth=0.3,
    #counts = True,
    #steps=100
    
).mark_area(opacity=0.5).encode(
    alt.X('ct:Q'),
    alt.Y('density:Q', stack='zero' ),
    alt.Color('Distribution_type:N')
)


chart.save(save_path+'density_plot_islet_vs_islet_control.html')


median_value = np.median(ROI_control_ct_distribution['ct'])

for i in range (0, len(file_name_list)):
    if ROI_control_ct_distribution['ct'][i]==median_value:
        print(file_name_list[i])

###############################################################################################################
'''

chart = alt.Chart(source).transform_fold(
    ['all_pairs', 'CCL19_CCR7'],
    as_=['Distribution Type', 'Attention Score']
).mark_bar(
    opacity=0.5,
    binSpacing=0
).encode(
    alt.X('Attention Score:Q', bin=alt.Bin(maxbins=100)),
    alt.Y('count()', stack=None),
    alt.Color('Distribution Type:N')
)

chart.save(save_path+'region_of_interest_filtered_combined_attention_distribution.html')
'''


############### Islet vs Acinar: find the 'ct' (and 'P'?) values for each of these ROI_control ########
ROI_control_Islet_Acinar_ct = dict() # column 3
ROI_control_ct_distribution = defaultdict(list)
file_name_list = []
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Acinar_Cells') or (out_histocat[i][1]=='Acinar_Cells' and out_histocat[i][2]=='Islet_Cells'):
            ROI_control_Islet_Acinar_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] == '1':
                ROI_control_ct_distribution['Distribution_type'].append('Islet vs Acinar cells in Control')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
            elif out_histocat[i][9] == 'NA': # or out_histocat[i][9] == '-1':
                continue
            else:
                ROI_control_ct_distribution['Distribution_type'].append('Islet vs Acinar cells in Control')
                ROI_control_ct_distribution['ct'].append(0)
            file_name_list.append(out_histocat[i][0])
#################### draw density plot of ct values for islet vs islet cells in control ###############


df = pd.DataFrame(ROI_control_ct_distribution)
chart =alt.Chart(df).transform_density(
    'ct',
    as_=['ct', 'density'],
    groupby=['Distribution_type'],
    #bandwidth=0.3,
    #counts = True,
    #steps=100
    
).mark_area(opacity=0.5).encode(
    alt.X('ct:Q'),
    alt.Y('density:Q', stack='zero' ),
    alt.Color('Distribution_type:N')
)

chart.save(save_path+'density_plot_islet_vs_acinar_control.html')
median_value = np.median(ROI_control_ct_distribution['ct'])

for i in range (0, len(file_name_list)):
    if ROI_control_ct_distribution['ct'][i]==median_value:
        print(file_name_list[i])


################################# Islet vs CD8+ cells in 4 states ###############################################################################
ROI_control_ct_distribution = defaultdict(list)

### in control ##
ROI_control_Islet_CD8_ct = dict() # column 3
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_control_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] == '1':
                ROI_control_ct_distribution['Distribution_type'].append('Control')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
            elif out_histocat[i][9] == 'NA':
                continue
            else:
                ROI_control_ct_distribution['Distribution_type'].append('Control')
                ROI_control_ct_distribution['ct'].append(0)


### in AAB+ ##
ROI_AAB_poz_Islet_CD8_ct = dict() # column 3
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_AAB_poz:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_AAB_poz_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] == '1':
                ROI_control_ct_distribution['Distribution_type'].append('AAB+')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
            elif out_histocat[i][9] == 'NA':
                continue
            else:
                ROI_control_ct_distribution['Distribution_type'].append('AAB+')
                ROI_control_ct_distribution['ct'].append(0)

### in Short T1DM ##
ROI_short_T1DM_Islet_CD8_ct = dict() # column 3
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_AAB_poz:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_short_T1DM_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] == '1':
                ROI_control_ct_distribution['Distribution_type'].append('T1DM<1')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
            elif out_histocat[i][9] == 'NA':
                continue
            else:
                ROI_control_ct_distribution['Distribution_type'].append('T1DM<1')
                ROI_control_ct_distribution['ct'].append(0)

### in Long T1DM ##
ROI_long_T1DM_Islet_CD8_ct = dict() # column 3
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_AAB_poz:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_long_T1DM_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] == '1':
                ROI_control_ct_distribution['Distribution_type'].append('T1DM>1')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
            elif out_histocat[i][9] == 'NA':
                continue
            else:
                ROI_control_ct_distribution['Distribution_type'].append('T1DM>1')
                ROI_control_ct_distribution['ct'].append(0)

####### make combined density plot ###########

df = pd.DataFrame(ROI_control_ct_distribution)
chart =alt.Chart(df).transform_density(
    'ct',
    as_=['ct', 'density'],
    groupby=['Distribution_type'],
    #bandwidth=0.3,
    #counts = True,
    #steps=100
    
).mark_area(opacity=0.5).encode(
    alt.X('ct:Q'),
    alt.Y('density:Q', stack='zero' ),
    alt.Color('Distribution_type:N')
)

chart.save(save_path+'density_plot_islet_vs_CD8poz_control_AAB_T1DM_short_long.html')
