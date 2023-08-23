import os
import pandas as pd
import copy
import csv
import numpy as np
import sys
import altair as alt
import altairThemes 
from collections import defaultdict

# register the custom theme under a chosen name
alt.themes.register("publishTheme", altairThemes.publishTheme)
# enable the newly registered theme
alt.themes.enable("publishTheme")

save_path = '/mnt/data1/fatema/'
###################### read cell labels, status, ROI filenames #######################################################

cell_label_all=[]
with open('/mnt/data1/fatema/cell_label_mgDF.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        cell_label_all.append(line[0])


cell_label=[]
with open('/mnt/data1/fatema/cell_label_islets_mgDF.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        cell_label.append(line[0])


x_coord=[]
with open('/mnt/data1/fatema/x_coord_mgDF.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        x_coord.append(line[0])

y_coord=[]
with open('/mnt/data1/fatema/y_coord_mgDF.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        y_coord.append(line[0])
        
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
        

######## find the list of regions for T1DM <1 or >1 status ###########
ROI_long_T1DM = dict()
ROI_short_T1DM = dict()
with open('/mnt/data1/gw/research/sc_integration/labels/labels_donor_t1dm_time.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        if 'T1DM>=1' == line[1]:
            ROI_long_T1DM[line[0]+'.tif']=''
        elif 'T1DM<1' == line[1]:
            ROI_short_T1DM[line[0]+'.tif']=''

###

###

## convert the dictionaries to list ##
ROI_control = list(ROI_control.keys())
ROI_AAB_poz = list(ROI_AAB_poz.keys())
ROI_long_T1DM = list(ROI_long_T1DM.keys())
ROI_short_T1DM = list(ROI_short_T1DM.keys())

########################### read histocat analysis result #############################################

out_histocat=[]
with open('/mnt/data1/fatema/out.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        out_histocat.append(line)

############### Islet vs Islet: find the 'ct' (and sigval?) values for each of these ROI_control ########

ROI_control_Islet_Islet_ct = dict() # column 3
box_plot_ct = defaultdict(list)
ROI_control_ct_distribution = defaultdict(list)
file_name_list = []
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Islet_Cells':
            
            if out_histocat[i][9] != 'NA': #= '1':
                ROI_control_ct_distribution['Distribution_type'].append('Islet vs Islet cells in Control')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
                box_plot_ct['y'].append(3.0)
                box_plot_ct['ct'].append(float(out_histocat[i][3]))
                file_name_list.append(out_histocat[i][0])
                ROI_control_Islet_Islet_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
                print('%s, %g'%(out_histocat[i][0], float(out_histocat[i][3])))
            
            else: #0, -1, NA
                continue

#################### draw density plot of ct values for islet vs islet cells in control ###############

df = pd.DataFrame(ROI_control_ct_distribution)
chart1 = alt.Chart(df).transform_density(
    'ct',
    as_=['ct', 'density'],
    groupby=['Distribution_type'],
    #bandwidth=0.3,
    #counts = True,
    #steps=100
    
).mark_area(opacity=0.7).encode(
    alt.X('ct:Q'),
    alt.Y('density:Q', stack='zero' ),
    alt.Color('Distribution_type:N'),
    tooltip=['ct']
)


#chart.save(save_path+'density_plot_islet_vs_islet_control_sigval1.html')
chart1.save(save_path+'density_plot_islet_vs_islet_control_sigval_any.html')

df = pd.DataFrame(box_plot_ct)
chart2 = alt.Chart(df).mark_boxplot(extent="min-max", opacity=0.7).encode(
    alt.X("ct:Q", stack='zero' ).scale(zero=False)
)
chart2.save(save_path+'box_plot_minmax_islet_vs_islet_control_sigval_any.html')

chart = alt.layer(chart1, chart2)
chart.save(save_path+'box_and_density_islet_vs_islet_control_sigval_any.html')



median_value = np.round(np.median(ROI_control_ct_distribution['ct']),2)
file_to_plot = ''
for i in range (0, len(file_name_list)):
    if np.round(ROI_control_ct_distribution['ct'][i],2)==median_value:
        print(file_name_list[i])
        file_to_plot = file_name_list[i]
        

##### draw the scatter plot ################
cell_meta_data = []
annospat_label = ['Acinar_Cells',  'Alpha_Cells', 'Beta_Cells', 'Delta_Cells',  'Gamma_Cells', 'Killer_T_Cells' ]
for cell_serial in range (0, len(file_name)):
    if file_name[cell_serial] == file_to_plot and cell_label_all[cell_serial] in annospat_label:
        cell_meta_data.append([float(x_coord[cell_serial]), float(y_coord[cell_serial]), cell_label_all[cell_serial] ])

data_list=dict()
data_list['X']=[]
data_list['Y']=[]   
data_list['label']=[] 
for i in range (0, len(cell_meta_data)):
    data_list['X'].append(cell_meta_data[i][0])
    data_list['Y'].append(cell_meta_data[i][1])
    data_list['label'].append(cell_meta_data[i][2])

set1 = altairThemes.get_colour_scheme("Set1", len(list(set(data_list['label']))))
#color=alt.Color('component_label:N', scale=alt.Scale(range=set1))
data_list_pd = pd.DataFrame(data_list)
chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
    alt.X('X', scale=alt.Scale(zero=False)),
    alt.Y('Y', scale=alt.Scale(zero=False)),
    color=alt.Color('label:N', scale=alt.Scale(range=set1)), #
    tooltip=['label'] #,'opacity'
)#.configure_legend(labelFontSize=6, symbolLimit=50)

#chart.save(save_path+'altair_plot_islet_islet_control_sigval1_median.html')
chart.save(save_path+'altair_plot_islet_islet_control_sigval_any_median.html')

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


############### Islet vs Acinar: find the 'ct' (and 'sigval'?) values for each of these ROI_control ########
ROI_control_Islet_Acinar_ct = dict() # column 3
ROI_control_ct_distribution = defaultdict(list)
file_name_list = []
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Acinar_Cells') or (out_histocat[i][1]=='Acinar_Cells' and out_histocat[i][2]=='Islet_Cells'):
            
            if out_histocat[i][9] != 'NA': # == '1':
                ROI_control_ct_distribution['Distribution_type'].append('Islet vs Acinar cells in Control')
                ROI_control_ct_distribution['ct'].append(float(out_histocat[i][3]))
                file_name_list.append(out_histocat[i][0])
                ROI_control_Islet_Acinar_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
                print('%s, %g'%(out_histocat[i][0], float(out_histocat[i][3])))
            else: #0, -1, NA
                continue

#################### draw density plot of ct values for islet vs acinar cells in control ###############

df = pd.DataFrame(ROI_control_ct_distribution)
chart1 =alt.Chart(df).transform_density(
    'ct',
    as_=['ct', 'density'],
    groupby=['Distribution_type']
    
).mark_area(opacity=0.7).encode(
    alt.X('ct:Q'),
    alt.Y('density:Q', stack='zero' ),
    alt.Color('Distribution_type:N'),
    tooltip=['ct']
)

#chart1.save(save_path+'density_plot_islet_vs_acinar_control_sigval1.html')
chart1.save(save_path+'density_plot_islet_vs_acinar_control_sigval_any.html')


chart2 = alt.Chart(df).mark_boxplot(extent="min-max", opacity=0.7).encode(
    alt.X("ct:Q", stack='zero').scale(zero=False)
)
chart2.save(save_path+'box_plot_minmax_islet_vs_acinar_control_sigval_any.html')

chart = alt.layer(chart1, chart2)
chart.save(save_path+'box_and_density_islet_vs_acinar_control_sigval_any.html')



median_value = np.round(np.median(ROI_control_ct_distribution['ct']),1)
file_to_plot = ''
for i in range (0, len(file_name_list)):
    if np.round(ROI_control_ct_distribution['ct'][i],1)==median_value:
        print(file_name_list[i])
        file_to_plot = file_name_list[i]
        

##### draw the scatter plot ################
cell_meta_data = []
annospat_label = ['Acinar_Cells',  'Alpha_Cells', 'Beta_Cells', 'Delta_Cells',  'Gamma_Cells', 'Killer_T_Cells' ]
for cell_serial in range (0, len(file_name)):
    if file_name[cell_serial] == file_to_plot and cell_label_all[cell_serial] in annospat_label:
        cell_meta_data.append([float(x_coord[cell_serial]), float(y_coord[cell_serial]), cell_label_all[cell_serial] ])

data_list=dict()
data_list['X']=[]
data_list['Y']=[]   
data_list['label']=[] 
for i in range (0, len(cell_meta_data)):
    data_list['X'].append(cell_meta_data[i][0])
    data_list['Y'].append(cell_meta_data[i][1])
    data_list['label'].append(cell_meta_data[i][2])

set1 = altairThemes.get_colour_scheme("Set1", len(list(set(data_list['label']))))
data_list_pd = pd.DataFrame(data_list)
chart = alt.Chart(data_list_pd).mark_point(filled=True, opacity = 1).encode(
    alt.X('X', scale=alt.Scale(zero=False)),
    alt.Y('Y', scale=alt.Scale(zero=False)),
    color=alt.Color('label:N', scale=alt.Scale(range=set1)),
    tooltip=['label'] #,'opacity'
)#.configure_legend(labelFontSize=6, symbolLimit=50)

#chart.save(save_path+'altair_plot_islet_acinar_control_sigval1_median.html')
chart.save(save_path+'altair_plot_islet_acinar_control_sigval_any_median.html')

################################# Islet vs CD8+ cells in 4 states ###############################################################################
ROI_ct_distribution = defaultdict(list)

### in control ##
ROI_control_ct_distribution = []
ROI_control_Islet_CD8_ct = dict() # column 3
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_control_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] !='NA': # == '1':
                ROI_ct_distribution['Distribution_type'].append('Control')
                ROI_ct_distribution['ct'].append(float(out_histocat[i][3]))
                ROI_control_ct_distribution.append(float(out_histocat[i][3]))
            else: 
                continue



### in AAB+ ##
ROI_AAB_poz_Islet_CD8_ct = dict() # column 3
ROI_AAB_poz_ct_distribution = []
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_AAB_poz:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_AAB_poz_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] !='NA': #  == '1':
                ROI_ct_distribution['Distribution_type'].append('AAB+')
                ROI_ct_distribution['ct'].append(float(out_histocat[i][3]))
                ROI_AAB_poz_ct_distribution.append(float(out_histocat[i][3]))
            else: 
                continue

### in Short T1DM ##
ROI_short_T1DM_Islet_CD8_ct = dict() # column 3
ROI_short_T1DM_ct_distribution = []
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_short_T1DM:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_short_T1DM_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] !='NA': #  == '1':
                ROI_ct_distribution['Distribution_type'].append('T1DM<1')
                ROI_ct_distribution['ct'].append(float(out_histocat[i][3]))
                ROI_short_T1DM_ct_distribution.append(float(out_histocat[i][3]))
            else: #out_histocat[i][9] == 'NA':
                continue

### in Long T1DM ##
ROI_long_T1DM_Islet_CD8_ct = dict() # column 3
ROI_long_T1DM_ct_distribution = []
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_long_T1DM:
        if  (out_histocat[i][1]== 'Islet_Cells' and out_histocat[i][2]== 'Killer_T_Cells') or (out_histocat[i][1]== 'Killer_T_Cells' and out_histocat[i][2]== 'Islet_Cells' ) :
            ROI_long_T1DM_Islet_CD8_ct[out_histocat[i][0]] = [out_histocat[i][3], out_histocat[i][9]]
            if out_histocat[i][9] !='NA': #  == '1':
                ROI_ct_distribution['Distribution_type'].append('T1DM>=1')
                ROI_ct_distribution['ct'].append(float(out_histocat[i][3]))
                ROI_long_T1DM_ct_distribution.append(float(out_histocat[i][3]))
            else: # out_histocat[i][9] == 'NA':
                continue

####### make combined density plot ###########
set1 = altairThemes.get_colour_scheme("Set1", len(list(set(ROI_ct_distribution['Distribution_type']))))
df = pd.DataFrame(ROI_ct_distribution)
chart1 =alt.Chart(df).transform_density(
    'ct',
    as_=['ct', 'density'],
    groupby=['Distribution_type'],    
).mark_area(opacity=0.7).encode(
    alt.X('ct:Q'),
    alt.Y('density:Q' ), #, stack='zero'
    alt.Color('Distribution_type:N', scale=alt.Scale(range=set1))
)
chart1.save(save_path+'density_plot_islet_vs_CD8poz_control_AAB_T1DM_short_long_sigval_any.html')
#chart1.save(save_path+'density_plot_islet_vs_CD8poz_control_AAB_T1DM_short_long_sigval1.html')

chart2 = alt.Chart(df).mark_boxplot(extent="min-max", opacity=0.7).encode(
    alt.X("ct:Q").scale(zero=False), 
    alt.Y('Distribution_type:N'),
    alt.Color('Distribution_type:N')
)
chart2.save(save_path+'box_plot_islet_vs_CD8poz_control_AAB_T1DM_short_long_sigval_any.html')

chart = alt.hconcat(chart1, chart2)
chart.save(save_path+'box_and_density_islet_vs_CD8poz_control_AAB_T1DM_short_long_sigval_any.html')

################### from scipy.stats import mannwhitneyu
# ROI_control_ct_distribution
# ROI_AAB_poz_ct_distribution
# ROI_short_T1DM_ct_distribution
# ROI_long_T1DM_ct_distribution

from scipy.stats import mannwhitneyu

dummy, p = mannwhitneyu(ROI_control_ct_distribution, ROI_AAB_poz_ct_distribution, method="exact")
print('mannwhitneyu p: control - AAB_poz = %g'%p)

dummy, p = mannwhitneyu(ROI_control_ct_distribution, ROI_short_T1DM_ct_distribution, method="exact")
print('mannwhitneyu p: control - short_T1DM = %g'%p)

dummy, p = mannwhitneyu(ROI_control_ct_distribution, ROI_long_T1DM_ct_distribution, method="exact")
print('mannwhitneyu p: control - long_T1DM = %g'%p)

dummy, p = mannwhitneyu(ROI_AAB_poz_ct_distribution, ROI_short_T1DM_ct_distribution, method="exact")
print('mannwhitneyu p: AAB_poz - short_T1DM =  %g'%p)

dummy, p = mannwhitneyu(ROI_AAB_poz_ct_distribution, ROI_long_T1DM_ct_distribution, method="exact")
print('mannwhitneyu p: AAB_poz - long_T1DM = %g'%p)

dummy, p = mannwhitneyu(ROI_short_T1DM_ct_distribution, ROI_long_T1DM_ct_distribution, method="exact")
print('mannwhitneyu p: short_T1DM - long_T1DM = %g'%p)

###############################################
from scipy.stats.mstats import kruskal

dummy, p = kruskal(ROI_control_ct_distribution, ROI_AAB_poz_ct_distribution)
print('kruskal p: control - AAB_poz = %g'%p)

dummy, p = kruskal(ROI_control_ct_distribution, ROI_short_T1DM_ct_distribution)
print('kruskal p: control - short_T1DM = %g'%p)

dummy, p = kruskal(ROI_control_ct_distribution, ROI_long_T1DM_ct_distribution)
print('kruskal p: control - long_T1DM = %g'%p)

dummy, p = kruskal(ROI_AAB_poz_ct_distribution, ROI_short_T1DM_ct_distribution)
print('kruskal p: AAB_poz - short_T1DM =  %g'%p)

dummy, p = kruskal(ROI_AAB_poz_ct_distribution, ROI_long_T1DM_ct_distribution)
print('kruskal p: AAB_poz - long_T1DM = %g'%p)

dummy, p = kruskal(ROI_short_T1DM_ct_distribution, ROI_long_T1DM_ct_distribution)
print('kruskal p: short_T1DM - long_T1DM = %g'%p)

kruskal(ROI_AAB_poz_ct_distribution, ROI_control_ct_distribution, ROI_short_T1DM_ct_distribution, ROI_long_T1DM_ct_distribution)
