"""
Created on Sat Jan 30 14:40:05 2021


@author: Carlos Caicedo-Montoya
"""

#Import libraries
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns



#Read data
my_data_piggy = pd.read_table('IGR_presence_absence_hygroscopicus.Rtab.txt',
                     sep='\t',
                     low_memory=False)

my_data_roary = pd.read_table('gene_presence_absence_hygroscopicus.Rtab.txt',
                     sep='\t',
                     low_memory=False)

# Change dataframe index
my_data_piggy.set_index("Gene", inplace=True)
my_data_roary.set_index("Gene", inplace=True)



# Pangenome frequency plot
colors = ["steelblue", "tan"]
clusters = []
clusters.append(np.asarray(my_data_roary.sum(axis=1), np.dtype(int)))
clusters.append(np.asarray(my_data_piggy.sum(axis=1), np.dtype(int)))
fig, ax = plt.subplots(figsize = (3.4, 3))

ax.hist(clusters, my_data_piggy.shape[1], histtype="bar", 
        color = colors, align= 'right')
#log = True)


ax.set_xlabel('Number of genomes', fontsize=10)
ax.set_ylabel('Number of clusters', fontsize=10)
ax.grid()
ax.legend(["Roary", "Piggy"])

#fig.savefig("roary-piggy-all.svg", dpi = 1200, format = 'svg')

my_data_roary["addition"] = my_data_roary.sum(axis=1)
frequency_roary = my_data_roary['addition'].value_counts()


my_data_piggy["addition"] = my_data_piggy.sum(axis=1)
frequency_piggy = my_data_piggy['addition'].value_counts()


