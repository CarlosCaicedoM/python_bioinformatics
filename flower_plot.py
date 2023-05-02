# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 19:52:32 2020

@author: Carlos Caicedo-Montoya
"""


#Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 

#Read data
my_data = pd.read_table('stats.txt',  sep='\t')
my_data.drop(["Genome no."], axis=1, inplace=True)
my_data.drop(list(my_data.columns[4:]), axis=1, inplace=True)


size = 0.3
unique_genes = my_data.iloc[:, 3]
accessory_genes = my_data.iloc[:, 2]
core_genes = my_data.iloc[:, 1]
core_genes = core_genes[0:2]
organisms = list(my_data.iloc[:,0])

cmap = plt.get_cmap("tab10") #### tab10, tab20, tab20b, tab20c


outer_colors = cmap(np.arange(len(unique_genes)))
inner_colors = cmap(np.arange(len(accessory_genes)))


fig, ax = plt.subplots(figsize = (8, 5))

ax.pie(unique_genes, radius=1, colors=outer_colors,
       labels = list(unique_genes), 
       labeldistance=0.79, frame = True, 
       wedgeprops=dict(width=size, edgecolor='w'))

ax.pie(accessory_genes, radius=1-size, colors=inner_colors,
       labels = list(accessory_genes), 
        labeldistance=0.59,  rotatelabels = 45,
       wedgeprops=dict(width=size, edgecolor='w'))

core_label = ["Core genome"]
core_label.append((str(core_genes[0])))

ax.pie(core_genes, radius=1-0.62, colors=['gold'], 
       labels = core_label, labeldistance=0.1)

ax.legend(organisms, loc = 'center right', fontsize = 10, 
          bbox_to_anchor=(1, 0, 0.5, 1))

ax.set(aspect="equal")


#fig.savefig("flower_plot", dpi = 1200, format = 'svg')



