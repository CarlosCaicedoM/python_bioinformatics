# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 21:43:25 2021

@author: Carlos Caicedo-Montoya
"""
import numpy as np
import pandas as pd

with open('BUSCO_concatenate.txt') as f:
    lines = f.readlines()

number_of_genomes = 121
Busco =[]
index = list(range(7, number_of_genomes*14, 14))
for i in index:
    A = lines[i]
    Busco.append(A)

#Get the complete list of genes
complete = [i.split(',F', 1)[0] for i in Busco]

##Get the fragmented list of genes
fragmented = [i.split(',', 4)[2] for i in Busco]
fragmented = [i.replace('F:', '') for i in fragmented]
fragmented = [i.replace('%', '') for i in fragmented]
fragmented = [float(i) for i in fragmented]
fragmented = np.asarray(fragmented)

##Get the missing list of genes
Missing = [i.split(',', 4)[3] for i in Busco]
Missing = [i.replace('M:', '') for i in Missing]
Missing = [i.replace('%', '') for i in Missing]
Missing = [float(i) for i in Missing]
Missing = np.asarray(Missing)

##Get the single list of genes
Single = [i.split('%', 3)[1] for i in complete]
Single =[i.replace('[S:', '') for i in Single]
Single = [float(i) for i in Single]
Single =  np.asarray(Single)

##Get the list of duplicated genes
Duplicated = [i.split('%', 3)[2] for i in complete]
Duplicated = [i.replace(',D:', '') for i in Duplicated]
Duplicated = [float(i) for i in Duplicated]
Duplicated = np.asarray (Duplicated)

#Complete
Cmplt = Single + Duplicated

#Create a dataframe with the results

df = pd.DataFrame({'Complete':Cmplt,
                   'Single':Single,
                   'Duplicated':Duplicated, 
                   'Fragmented':fragmented,
                   'Missing':Missing})


#plot the results
import seaborn as sns
import matplotlib.pyplot as plt
 
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(7, 7))
with sns.axes_style("white"):
    sns.histplot( df["Single"] , color="skyblue", ax=ax1)
    sns.histplot( df["Duplicated"] , color="olive", ax=ax2)
    sns.histplot( df["Fragmented"] , color="gold", ax=ax3)
    sns.histplot( df["Missing"] , color="teal", ax=ax4)
ax1.set_ylabel("Number of genomes")
ax2.set_ylabel("")
ax3.set_ylabel("Number of genomes")
ax4.set_ylabel("")
ax1.set_xlabel("")
ax2.set_xlabel("")
ax3.set_xlabel("%BUSCOs")
ax4.set_xlabel("%BUSCOs")

ax1.set_title("Complete and single-copy")
ax2.set_title("Complete and duplicated")
ax3.set_title("Fragmented")
ax4.set_title("Missing")

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()


f.savefig("FigureS1", dpi = 1200, format = 'pdf')
