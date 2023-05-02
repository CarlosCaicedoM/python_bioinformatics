# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 23:30:01 2021

@author: User
"""
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

my_data = pd.read_table("GID.txt", sep='\t', low_memory=False)

GC = my_data.iloc[:,4]
GC = np.asarray(GC)

print(np.mean(GC))
print(np.std(GC))


fig1, ax = plt.subplots()
with sns.axes_style("whitegrid"):
    sns.histplot(GC, color="darkblue", ax=ax, binwidth=.5)
ax.grid()
ax.set_xlabel('%GC')
ax.set_ylabel('Number of genomes') 


#Genome size and number of proteins correlation
y1 = my_data.iloc[:,2]
y1 = np.asarray(y1)   #Sequence length
x1 = my_data.iloc[:,3]
x1 = np.asarray(x1)   #Number of proteins


slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)


x = np.linspace(5000, 10000, 1000) 
my_line = slope*x+intercept


#Genome size and correlation with GC contentc
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(y1, GC)

slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(GC, x1)


#plot the results
f, ax = plt.subplots()
ax.plot(x1, y1, 'oy')
ax.plot(x, my_line, '-b')
ax.text(6000, 10.5, r'$R^2 = {0:3.2f}$'.format(r_value**2))
ax.text(5900, 10,
         r'$y={0:3.4f}x + {1:3.3f}$'.format(slope, intercept))
ax.grid()
ax.set_xlabel("Number of proteins")
ax.set_ylabel("Genome size (Mb)")

#3D PLOT
fig1 = plt.figure()
ax1 = plt.axes(projection="3d")
ax1.scatter3D(GC, x1, y1, c=y1, cmap='RdYlBu')
ax1.view_init(20, 60)
ax1.set_xlabel('%G+C')
ax1.set_ylabel('Number of proteins')
ax1.set_zlabel('Genome size (Mb)')

fig2 = plt.figure()
ax2 = plt.axes(projection="3d")
ax2.scatter3D(GC, x1, y1, c=y1, cmap='viridis')
ax2.set_xlabel('%G+C')
ax2.set_ylabel('Number of proteins')
ax2.set_zlabel('Genome size (Mb)')

fig2.savefig("FigureS2", dpi = 1200, format = 'pdf')

