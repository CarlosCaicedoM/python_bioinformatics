# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 08:58:11 2020

@author: Carlos Caicedo-Montoya

Tutorial Heatmap in Python
"""

from bioinfokit import analys, visuz

df = analys.get_data('hmap').data


# set gene names as index
df = df.set_index(df.columns[0])

# heatmap with hierarchical clustering 
visuz.gene_exp.hmap(df=df, dim=(3, 6), tickfont=(6, 4))

# heatmap without hierarchical clustering 
visuz.gene_exp.hmap(df=df, dim=(3, 6), tickfont=(6, 4), rowclus=False, show=True)


# colormaps are available at  https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
# default is seismic 
# here I use red-yellow-green: RdYlGn
visuz.gene_exp.hmap(df=df, cmap='RdYlGn', dim=(3, 6), tickfont=(6, 4), show=True)

# Z-score can be used to standardize value with mean 0 and var 1
# default Z-score is set to None and it applies to only heatmap with cluster
# here I standardize column with Z-score
visuz.gene_exp.hmap(df=df, zscore=0, dim=(3, 6), tickfont=(6, 4), show=True)


from bioinfokit import analys, visuz
import pandas as pd
import os


my_data =pd.read_table("Streptomyces_Ajit.txt")


# set gene names as index
my_data = my_data.set_index(my_data.columns[0])

# heatmap with hierarchical clustering 
visuz.gene_exp.hmap(df=my_data, dim=(3, 6), tickfont=(6, 4))


#Roary Pangenome plots

#plottinh imports
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')

#other imports
import os
import pandas as pd
import numpy as np
from Bio import Phylo


#parSNP tree
#Any other valid newick file is fine, if the tip labels is
#the same as in the gene_presence_absence matrix from roary.

t = Phylo.read('Accessory_Binary_Genes.nwk', 'newick')

# Max distance to create better plots
mdist = max([t.distance(t.root, x) for x in t.get_terminals()])


# Load roary
roary = pd.read_table('gene_presence_absence.csv',
                     sep=',',
                     low_memory=False)


# Set index (group name)
roary.set_index('Gene', inplace=True)

# Drop the other info columns
roary.drop(list(roary.columns[:13]), axis=1, inplace=True)

# Transform it in a presence/absence matrix (1/0)
roary.replace('.{2,100}', 1, regex=True, inplace=True)
roary.replace(np.nan, 0, regex=True, inplace=True)


# Sort the matrix by the sum of strains presence
idx = roary.sum(axis=1).sort_values(ascending=False).index
roary_sorted = roary.loc[idx]

# Pangenome frequency plot
plt.figure(figsize=(7, 5))

plt.hist(roary.sum(axis=1), roary.shape[1],
         histtype="stepfilled", alpha=.7)

plt.xlabel('Number of genomes')
plt.ylabel('Number of genes')

sns.despine(left=True,
            bottom=True)


# Sort the matrix according to tip labels in the tree
roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]


# PLot presence/absence matrix against the tree
with sns.axes_style('whitegrid'):
    fig = plt.figure(figsize=(17, 10))

    ax1=plt.subplot2grid((1,40), (0, 10), colspan=30)
    a=ax1.imshow(roary_sorted.T, cmap=plt.cm.Blues,
               vmin=0, vmax=1,
               aspect='auto',
               interpolation='none',
                )
    ax1.set_yticks([])
    ax1.set_xticks([])
    ax1.axis('off')

    ax = fig.add_subplot(1,2,1)
    ax=plt.subplot2grid((1,40), (0, 0), colspan=10, facecolor='white')

    fig.subplots_adjust(wspace=0, hspace=0)

    ax1.set_title('Roary matrix\n(%d gene clusters)'%roary.shape[0])
    
    Phylo.draw(t, axes=ax, 
               show_confidence=False,
               label_func=lambda x: None,
               xticks=([],), yticks=([],),
               ylabel=('',), xlabel=('',),
               xlim=(-0.01,mdist+0.01),
               axis=('off',),
               title=('parSNP tree\n(%d strains)'%roary.shape[1],), 
              )
