# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 22:26:22 2021

@author: Carlos Caicedo-Montoya
"""

#Import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import Phylo
import numpy as np
import math

#Read_data
my_matrix = pd.read_csv('pm.pfam.table2_file_Streptomyces_pangenomeR.RData.csv', 
                        low_memory=False)
my_matrix.set_index('Genome', inplace=True)


my_matrix = my_matrix.T

# Transform it in a presence/absence matrix (1/0)
max_value= max(list(my_matrix.max()))
presence_absence_matrix = my_matrix.replace(list(range(1,max_value+1)), 1)
max(list(presence_absence_matrix.max()))



genomes = presence_absence_matrix.columns.tolist()

genome_size = my_matrix.sum(axis=0) 
genome_size_per_domains = presence_absence_matrix.sum(axis=0) 



core_superior = math.ceil(len(presence_absence_matrix.T))
core_inferior = math.ceil(len(presence_absence_matrix.T)*0.95)
shell_inferior = math.ceil(len(presence_absence_matrix.T)*0.15)

genome_core  = {}
genome_shell  = {}
genome_cloud  = {}

addition = presence_absence_matrix.sum(axis=1) 
presence_absence_matrix['addition'] = addition



soft_core = presence_absence_matrix[presence_absence_matrix["addition"] >= core_inferior]
shell = presence_absence_matrix[(presence_absence_matrix["addition"] >= shell_inferior) & (presence_absence_matrix["addition"] < core_inferior)]
cloud = presence_absence_matrix[presence_absence_matrix["addition"] < shell_inferior]





soft_core_clusters = pd.DataFrame(list(soft_core.index))
shell_clusters = pd.DataFrame(list(shell.index))
cloud_clusters = pd.DataFrame(list(cloud.index))


soft_core_clusters.to_csv ('soft_core_clusters.csv', index = False, header=False)
shell_clusters.to_csv ('shell_clusters.csv', index = False, header=False)
cloud_clusters.to_csv ('cloud_clusters.csv', index = False, header=False)



