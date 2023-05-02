# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 13:43:03 2020


@author: Carlos Caicedo-Montoya
"""

#Import libraries
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



#Read data
my_data = pd.read_table('matrix_streptomyces_annotation.txt',
                     sep='\t',
                     low_memory=False)



# Drop the other info columns
my_data.drop(["Cluster", "GI", "Annotation"], axis=1, inplace=True)
my_data.drop(list(my_data.columns[-2:]), axis = 1, inplace = True)

# Pangenome frequency plot
fig, ax = plt.subplots(figsize = (3.4, 3))
ax.hist(my_data.sum(axis=1), my_data.shape[1], histtype="stepfilled",
         alpha=1, color = 'darkred')

ax.set_xlabel('Number of genomes', fontsize=10)
ax.set_ylabel('Number of genes', fontsize=10)


# Pie chart, where the slices will be ordered and plotted counter-clockwise:
core = my_data[my_data.sum(axis=1) == my_data.shape[1]].shape[0]
unique = my_data [my_data.sum(axis = 1) ==1].shape[0]
accessory = my_data.shape[0] - core - unique


     
labels = ['Core \n({0:0.0f})'.format(core), 
          'Accessory \n({0:0.0f})'.format(accessory), 
          'Unique \n ({0:0.0f})'.format(unique)]

sizes = [core, accessory, unique]
explode = (0.1, 0.051, 0.05)  

fig1, ax1 = plt.subplots(figsize = (3.4, 3))
ax1.pie(sizes, explode=explode, labels=labels, shadow=True, 
        startangle = 90, autopct='%0.1f%%', 
        colors = ["cornflowerblue", "lightsalmon", "turquoise"])
        
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

fig1.savefig("pie_chart_streptomyces", dpi = 1200, format = 'svg')