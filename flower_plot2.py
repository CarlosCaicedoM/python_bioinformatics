# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 10:48:36 2020


@author: Carlos Caicedp-Montoya

"""


#Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 

#Read data
my_data = pd.read_table('stats_streptomyces.txt',  sep='\t')
my_data.drop(["Genome no."], axis=1, inplace=True)
my_data.drop(list(my_data.columns[4:]), axis=1, inplace=True)


size = 0.3
unique_genes = my_data.iloc[:, 3]
accessory_genes = my_data.iloc[:, 2]
core_genes = my_data.iloc[:, 1]
organisms = list(my_data.iloc[:,0])


petals = np.ones(len(unique_genes))


fig, ax = plt.subplots(figsize = (8, 8))
ax.pie(petals, radius=1.2, colors=["turquoise"],
       labels = organisms, 
       labeldistance=1.05, frame = True, 
       rotatelabels = 45, textprops=dict(fontsize=14),
       wedgeprops=dict(width=size, edgecolor='w'))

ax.pie(petals, radius=0.9, colors=["lightsalmon"],
       labels = list(unique_genes), 
        labeldistance=1,  rotatelabels = 45,
       wedgeprops=dict(width=size, edgecolor='w'), textprops=dict(fontsize=14))

ax.pie(core_genes, radius=0.57, colors=['cornflowerblue'], 
       labels = list(accessory_genes), textprops=dict(fontsize=12),
       rotatelabels = 45, labeldistance=1.051,
       wedgeprops=dict(width=1))

ax.text(-0.23, 0, 'Core genome ', fontsize = 16)
ax.text(-0.1, -0.2, '{0:0.0f}'.format(core_genes[0]), fontsize = 16)
ax.set(aspect="equal")


#fig.savefig("flower_plot_Streptomyces2", dpi = 1200, format = 'svg')

"""
'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG',
 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu',
 'BuPu_r', 'CMRmap', 'CMRmap_r', 
 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r',
 'Greens', 'Greens_r', 'Greys', 'Greys_r', '
 OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 
 'PRGn', 'PRGn_r', 'Paired', 'Paired_r',
 'Pastel1', 'Pastel1_r', 'Pastel2', 
 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r',
 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd',
 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 
 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 
 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 
 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 
 'Spectral', 'Spectral_r', 'Wistia',
 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 
 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r',
 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 
 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 
 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 
 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix',
 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 
 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat',
 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow',
 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg',
 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r',
 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 
 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma',
 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 
 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 
 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r',
 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer',
 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 
 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 
 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r',
 'twilight_shifted', 'twilight_shifted_r', 'viridis', 
 'viridis_r', 'winter', 'winter_r'
"""