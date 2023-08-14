import os
#import glob
import pandas as pd
#import shutil
import copy
import csv
import numpy as np
import sys

###########################################################################################

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

# find the regions for 'control' or 'GAD-' status
ROI_control = dict()
for cell_serial in range (0, len(file_name)):
    if status_list[cell_serial] == 'GAD-':
        ROI_control[file_name[cell_serial]] = ''

print(ROI_control.keys())
print('count of ROI with control state is %d'%len(ROI_control.keys()))
ROI_control = list(ROI_control.keys())

#######################################

out_histocat=[]
with open('/mnt/data1/fatema/out.csv') as file:
    csv_file = csv.reader(file, delimiter=",")
    for line in csv_file:
        out_histocat.append(line)

# find the 'ct' (and 'P'?) values for each of these ROI_control
ROI_control_Islet_Islet_ct = dict() # column 3
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  == 'Islet_Cells' and == 'Islet_Cells':
        ROI_control_Islet_Islet_ct[out_histocat[i][0]] = out_histocat[i][3]

ROI_control_Islet_Acinar_ct = dict() # column 3
for i in range (1, len(out_histocat)):
    if out_histocat[i][0] in ROI_control:
        if  == 'Islet_Cells' and == 'Islet_Cells':
        ROI_control_Islet_Acinar_ct[out_histocat[i][0]] = out_histocat[i][3]

df = pd.DataFrame(ROI_control_Islet_Islet_ct)
df.to_csv('mnt/data1/fatema/ROI_control_Islet_Islet_ct_mgDF.csv', index=False, header=False)

df = pd.DataFrame(ROI_control_Islet_Acinar_ct)
df.to_csv('mnt/data1/fatema/ROI_control_Islet_Acinar_ct_mgDF.csv', index=False, header=False)
##############################################
# do the plotting 
    

