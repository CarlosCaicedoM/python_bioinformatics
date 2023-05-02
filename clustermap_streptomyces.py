# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 13:49:02 2021

@author: User
"""
    
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np

my_matrix = pd.read_table("ANI_matrix.csv", sep = ',', 
                        low_memory = False)

my_matrix.set_index('reference', inplace=True)


my_colors = pd.read_table("colors.csv", sep = ',', 
                        low_memory = False)
my_colors.set_index('strain', inplace=True)


#Plot
sns_plot = sns.clustermap(my_matrix, linewidths=0.0, 
            cmap = 'viridis',
            dendrogram_ratio=(0.1, 0.2), 
            cbar_pos=(0.5, 0, .4, 0.03), 
            row_colors=my_colors, 
            col_colors=my_colors, 
            cbar_kws={"orientation": "horizontal"},
            tree_kws = {"linewidths":1.7},
            annot_kws = {'fontsize':8}, 
            figsize=(40, 40))
        
sns_plot.ax_col_dendrogram.remove()

#metric = "correlation"
#method="single"
#row_cluster = False---> Eliminate horizontal clustering
#col_cluster=False----Eliminate vertical clustering
#dendrogram_ratio=0.2; default

#figsize=(20, 20))

sns_plot.savefig("output4", dpi = 1200, format = 'svg')
#tree_kws = {"facecolors":list('rgb'), 
                        #"linewidths":1.6},