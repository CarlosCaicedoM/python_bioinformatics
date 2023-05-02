# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:19:58 2021

@author: Carlos Caicedo-Montoya
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import Phylo
import numpy as np
import math

#Read_data
my_matrix = pd.read_csv('pm.pfam.table2_file_Streptomyces_pangenomeR.RData.csv',
                        low_memory=False)

# Set index (group name)
my_matrix.set_index('Genome', inplace=True)

#Transpose pfam matrix
#As produced by Micropan, the matrix must be transposed to use similar
#codes to analyze the matrix in a similar way as in the case of BPGA and Roary

my_matrix = my_matrix.T

# Transform it in a presence/absence matrix (1/0)
max_value= max(list(my_matrix.max()))
presence_absence_matrix = my_matrix.replace(list(range(1,max_value+1)), 1)
max(list(presence_absence_matrix.max()))



## Pangenome frequency plot
core_superior = math.ceil(len(presence_absence_matrix.T))
core_inferior = math.ceil(len(presence_absence_matrix.T)*0.95)
shell_inferior = math.ceil(len(presence_absence_matrix.T)*0.15)

fig, ax = plt.subplots(figsize = (7, 7))
ax.hist(presence_absence_matrix.sum(axis=1),
        presence_absence_matrix.shape[1], histtype="stepfilled",
         alpha=1, color = 'crimson')

ax.set_xlabel('Number of genomes', fontsize=10)
ax.set_ylabel('Number of clusters', fontsize=10)
ax.grid()
ax.annotate('',
            xy=(19, 500), xycoords='data',
            xytext=(25, 40), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.01))
ax.annotate('',
            xy=(60, 100), xycoords='data',
            xytext=(0, 40), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.01))
ax.annotate('',
            xy=(120, 900), xycoords='data',
            xytext=(-20, 40), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.01))


fig.savefig("hist_Micropan", dpi = 1200, format = 'svg')


# Pie chart, where the slices will be ordered and plotted counter-clockwise:
softcore = presence_absence_matrix[presence_absence_matrix.sum(axis=1) >= core_inferior].shape[0]
shell = presence_absence_matrix[(presence_absence_matrix.sum(axis=1) >= shell_inferior) & (presence_absence_matrix.sum(axis=1) < core_inferior)].shape[0]
cloud = presence_absence_matrix.shape[0] - softcore - shell
     
labels = ['Shell \n({0:0.1f})'.format(shell),
          'Softcore \n({0:0.1f})'.format(softcore),  
          'Cloud \n ({0:0.1f})'.format(cloud)]

sizes = [shell, softcore, cloud]
explode = (0.1, 0.1, 0.0)  

fig1, ax1 = plt.subplots(figsize = (7, 7))
ax1.pie(sizes, explode=explode, labels=labels, shadow=True, 
        startangle = 220, autopct='%0.1f%%', 
        colors = ["gold", "darkseagreen", "skyblue"])
        
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

fig1.savefig("pie_chart_Micropan", dpi = 1200, format = 'svg')
