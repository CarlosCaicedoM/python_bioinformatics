# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 23:33:14 2021

@author: Carlos Caicedo

Establish the nmber of core, shell and cloud genes in each genome
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


core_superior = math.ceil(len(my_matrix.T))
core_inferior = math.ceil(len(my_matrix.T)*0.95)
shell_inferior = math.ceil(len(my_matrix.T)*0.15)

genomes = my_matrix.columns.tolist()

genome_size = my_matrix.sum(axis=0) 

genome_core  = {}
genome_shell  = {}
genome_cloud  = {}


addition = my_matrix.sum(axis=1) 
my_matrix['addition'] = addition


soft_core = my_matrix[my_matrix["addition"] >= core_inferior]
shell = my_matrix[(my_matrix["addition"] >= shell_inferior) & (my_matrix["addition"] < core_inferior)]
cloud = my_matrix[my_matrix["addition"] < shell_inferior]

for i in genomes:
    genome_core[i] = (sum(soft_core[i])) 
for i in genomes:
    genome_shell[i] = (sum(shell[i])) 
for i in genomes:
    genome_cloud[i] = (sum(cloud[i])) 



my_results = pd.DataFrame({'genome_core':pd.Series(genome_core),
              'genome_shell':pd.Series(genome_shell),
              'genome_cloud':pd.Series(genome_cloud)})


genome_size2 = my_results.sum(axis=1) 

genome_size2.equals(genome_size)

my_results.to_csv ('stats_core_shell_cloud.csv', index = True, header=True)




