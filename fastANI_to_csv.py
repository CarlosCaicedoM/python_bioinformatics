# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 14:04:38 2021

@author: User
"""
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt

my_data = pd.read_table("fastANI.txt", sep = '\t', 
                        low_memory = False)

my_matrix = my_data.set_index(['query', 'reference']).unstack('reference')

my_matrix.to_csv ('ANI_matrix.csv',
                    index = True, header=True)
