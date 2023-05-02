# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 19:21:23 2021

@author: CArlos Caicedo-Montoya
"""

#Import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import Phylo
import numpy as np
import math

#Read_data
my_matrix = pd.read_csv('Clusters_pfam_names.txt', 
                        sep='\t',
                        low_memory=False)

my_matrix.drop_duplicates(subset = "cluster", inplace = True)
my_matrix.to_csv ('Representative_sequences', index = False, header=True)


