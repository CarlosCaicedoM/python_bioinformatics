# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 23:33:14 2021

@author: Carlos Caicedo

Establish the nmber of core, accessory and unique genes in each genome
from a pangenome matrix

"""
#Import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import Phylo
import numpy as np
import math

#Read_data
my_matrix = pd.read_csv('roary_presence_absence.csv', low_memory=False)

# Set index (group name)
my_matrix.set_index('Gene', inplace=True)

# Drop the other info columns
my_matrix.drop(list(my_matrix.columns[:13]), axis=1, inplace=True)
#my_matrix.drop('addition', axis=1, inplace=True)

genomes = my_matrix.columns.tolist()

genome_size = my_matrix.sum(axis=0) 

genome_core  = {}
genome_accessory  = {}
genome_unique  = {}


addition = my_matrix.sum(axis=1) 
my_matrix['addition'] = addition


core = my_matrix[my_matrix["addition"] == max(addition)]
accessory = my_matrix[(my_matrix["addition"] > min(addition)) & (my_matrix["addition"] < max(addition))]
unique = my_matrix[my_matrix["addition"] == min(addition)]

for i in genomes:
    genome_core[i] = (sum(core[i])) 
for i in genomes:
    genome_accessory[i] = (sum(accessory[i])) 
for i in genomes:
    genome_unique[i] = (sum(unique[i])) 



my_results = pd.DataFrame({'genome_core':pd.Series(genome_core),
              'genome_shell':pd.Series(genome_accessory),
              'genome_cloud':pd.Series(genome_unique)})


genome_size2 = my_results.sum(axis=1) 

genome_size2.equals(genome_size)

my_results.to_csv ('stats_core_accessory_unique.csv', index = True, header=True)




