# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 15:30:04 2021

@author: Carlos Caicedo-Montoya
"""

#Pangenome plots

#plottinh imports
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns


#other imports
import os
import pandas as pd
import numpy as np



# Load BCG results
BCG = pd.read_table('BCG_overall.txt',
                     sep='\t', low_memory = False)


fig, ((ax1, ax2, ax3, ax4, ax5)) = plt.subplots(nrows=5, sharex = True, figsize=(15, 20))
sns.set_theme(style="whitegrid")
sns.barplot(x="BCG", y="Presence", data=BCG, ax=ax1)
sns.barplot(x="BCG", y="Core", data=BCG, ax=ax2)
sns.barplot(x="BCG", y="Res", data=BCG, ax=ax3)
sns.barplot(x="BCG", y="BGC_Phyl_Res_Dupl", data=BCG, ax=ax5)
sns.barplot(x="BCG", y="BGC_Phyl_Res", data=BCG, ax=ax4)

ax5.tick_params(top=False, bottom=True, which = 'both',
                  labeltop=False, labelbottom=True)


plt.setp(ax1.get_xticklabels(), rotation=90, ha="center",
         rotation_mode="default")



fig.savefig("test3.svg", dpi = 1200, format = "svg")