# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 23:04:58 2021

@author: User
"""
import pandas as pd
from Bio import Phylo
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns

#Read_data_Roary
my_matrix_roary= pd.read_csv('roary_presence_absence.csv', low_memory=False)

# Set index (group name)
my_matrix_roary.set_index('Gene', inplace=True)

# Drop the other info columns
my_matrix_roary.drop(list(my_matrix_roary.columns[:13]), axis=1, inplace=True)
#my_matrix.drop('addition', axis=1, inplace=True)

#Read_data_Micropan
my_matrix_micropan = pd.read_csv('pm.pfam.table2_file_Streptomyces_pangenomeR.RData.csv',
                        low_memory=False)

# Set index (group name)
my_matrix_micropan.set_index('Genome', inplace=True)

#Transpose pfam matrix
#As produced by Micropan, the matrix must be transposed to use similar
#codes to analyze the matrix in a similar way as in the case of BPGA and Roary

my_matrix_micropan = my_matrix_micropan.T

# Transform it in a presence/absence matrix (1/0)
max_value= max(list(my_matrix_micropan.max()))
my_matrix_micropan = my_matrix_micropan.replace(list(range(1,max_value+1)), 1)

## Pangenome frequency plot
core_superior = math.ceil(len(my_matrix_roary.T))
core_inferior = math.ceil(len(my_matrix_roary.T)*0.95)
shell_inferior = math.ceil(len(my_matrix_roary.T)*0.15)

#Pangenome size Roary
softcore_roary = my_matrix_roary[my_matrix_roary.sum(axis=1) >= core_inferior].shape[0]
shell_roary = my_matrix_roary[(my_matrix_roary.sum(axis=1) >= shell_inferior) & (my_matrix_roary.sum(axis=1) < core_inferior)].shape[0]
cloud_roary = my_matrix_roary.shape[0] - softcore_roary- shell_roary

#PAngenome size micropan
softcore_micropan = my_matrix_micropan[my_matrix_micropan.sum(axis=1) >= core_inferior].shape[0]
shell_micropan = my_matrix_micropan[(my_matrix_micropan.sum(axis=1) >= shell_inferior) & (my_matrix_micropan.sum(axis=1) < core_inferior)].shape[0]
cloud_micropan = my_matrix_micropan.shape[0] - softcore_micropan - shell_micropan

#plot
labels = ["soft-core", "shell", "cloud"]
roary =[softcore_roary, shell_roary, cloud_roary]
micropan = [softcore_micropan, shell_micropan, cloud_micropan] 
width = 0.4

fig, ax = plt.subplots(figsize=(7,7))
ax.bar(labels, roary, width, label='Roary', color = "steelblue")
ax.bar(labels, micropan, width, bottom=roary, label='Micropan', color ="orange")
ax.set_ylabel("Number of clusters")
ax.grid()
ax.legend()
ax.annotate('{}'.format(softcore_roary),
                    xy=(-0.15, 3000),
                    xytext=(-0.5, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', 
                    arrowprops=dict(facecolor='steelblue', shrink=0.01))
ax.annotate('{}'.format(softcore_micropan),
                    xy=(0.15, 3000),
                    xytext=(-0.5, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', 
                    arrowprops=dict(facecolor='orange', shrink=0.01))
ax.annotate('{}'.format(shell_roary),
                    xy=(0.85, 8000),
                    xytext=(-0.5, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', 
                    arrowprops=dict(facecolor='steelblue', shrink=0.01))
ax.annotate('{}'.format(shell_micropan),
                    xy=(1.15, 8000),
                    xytext=(-0.5, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', 
                    arrowprops=dict(facecolor='orange', shrink=0.01))
ax.annotate('{}'.format(cloud_roary),
                    xy=(1.85, 143000),
                    xytext=(-0.5, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', 
                    arrowprops=dict(facecolor='steelblue', shrink=0.01))
ax.annotate('{}'.format(cloud_micropan),
                    xy=(2.15, 143000),
                    xytext=(-0.5, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', 
                    arrowprops=dict(facecolor='orange', shrink=0.01))

fig.savefig("pan_sizes", dpi=1200)

