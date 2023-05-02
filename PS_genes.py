# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 17:46:25 2020

@author: Carlos Caicedo-Montoya
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#READ DATA
my_data = pd.read_table('COG_summary.txt',
                     sep='\t',
                     low_memory=False)


labels = my_data.iloc[:, 0]
labels = list(labels)

L06 = my_data.loc[:,"Abundance_L06"]
L06 = np.asarray(L06)
L06 = L06*100

K155 = my_data.loc[:,"Abundance_K155"]
K155 = np.asarray(K155)

NF3 = my_data.loc[:,"Abundance_NF3"]
NF3 = np.asarray(NF3)
NF3 = NF3*100

TFC3 = my_data.loc[:,"Abundance_TFC3"]
TFC3 = np.asarray(TFC3)
TFC3 = TFC3*100
x = np.arange(len(labels))  # the label locations
width = 0.2  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width, L06, width, label='L06')
rects2 = ax.bar(x, K155, width, label='K155')
rects3 = ax.bar(x + width, NF3, width, label='NF3')
#rects4 = ax.bar(x + width, TFC3, width, label='TFC3')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Percentage of total PS genes')
#ax.set_title('Scores by group and gender')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
#ax.grid()

