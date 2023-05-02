# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 23:58:24 2020

@author: Carlos Caicedo-Montoya
Heat map for gene copy number

"""

#Import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import Phylo
import numpy as np

#TREE
#newick file
#the tip labels should be  the same as in the gene_presence_absence matrix.

#my_tree = Phylo.read('MEGA-Alignment_panTree.nwk', 'newick')
my_tree = Phylo.read('MAFFT_core_phylogeny_2644_strain_names.nwk', 'newick')

# Load matrix
#archives
#heat_map_oxidative_stress
#heat_map_illumination2

my_matrix = pd.read_table('heavy_metal2.txt',
                     sep='\t',
                     low_memory=False)


# Set index (group name)
my_matrix.set_index('strain', inplace=True) #inplace parameter will change the dataframe without assignment


my_matrix = my_matrix.transpose()

# Sort the matrix according to tip labels in the tree
my_matrix_sorted = my_matrix[[x.name for x in my_tree.get_terminals()]]


fig1 = plt.figure(figsize=(15, 6))
ax1 = plt.subplot2grid((1, 3), (0, 1), colspan=2, fig = fig1)
sns.heatmap(my_matrix_sorted.T, linewidths=0.5, ax=ax1, 
            cmap = plt.cm.RdYlGn,
            cbar_kws = dict(ticks=np.arange(0, 15)), 
            yticklabels=False)
ax1.tick_params(top=True, bottom=False, which = 'both',
                   labeltop=True, labelbottom=False)


plt.setp(ax1.get_xticklabels(), rotation=90, ha="center",
         rotation_mode="default")


ax2=plt.subplot2grid((1,3), (0, 0), colspan=1, facecolor='white')
#f.subplots_adjust(wspace=0, hspace=0)
Phylo.draw(my_tree, axes=ax2, show_confidence=False, 
               xticks=([],), yticks=([],),
               ylabel=('',), xlabel=('',),
               xlim=(-0.1, 0.14),
               axis=('off',))

fig1.savefig("heavy_metal2", dpi = 1200, format = "svg")
   

#Colors

##Illumination  #coolwarm
##Oxidative stress #Spectral
##heavy meta #RdYlGn
##cell_wall, hyperosmotics  # RdYlBu


