# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 01:47:45 2020

@author: User
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import Phylo
import numpy as np


# Load matrix
my_matrix = pd.read_table('Actinoplanes_f0_0taxa_CDS_algOMCL_e0_Avg_identity_names.txt',
                     sep='\t',
                     low_memory=False)

my_matrix.set_index('genomes', inplace=True) #inplace parameter will change the dataframe without assignment



sns_plot = sns.clustermap(my_matrix, linewidths=0.5, 
            cmap = plt.cm.Spectral)

#sns_plot.savefig("output.png", dpi = 1200)
