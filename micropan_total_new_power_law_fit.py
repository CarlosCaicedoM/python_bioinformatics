# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 14:02:56 2021

@author: Usuario
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 13:11:38 2021
@author: Carlos Caicedo-Montoya
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

#Read data Pangenome
my_total_genes = pd.read_table("pangenome_micropan.csv",
                                   sep=",", 
                                   low_memory=False, header=None)
my_total_genes = my_total_genes.T

total_genes_melted = my_total_genes.melt()


genomes = total_genes_melted.iloc[:,0]
genomes = np.asarray(genomes)
genomes = genomes + 1
pangenome = total_genes_melted.iloc[:,1]


#Fit power law curve
#PANGENOME

def pangenome_fit(x, gamma, G0):
    return G0*(x**gamma)
parameters_pan, parameters_cov_pan = curve_fit(pangenome_fit,
                                           genomes,
                                           pangenome, p0=[1, 1])

#Generate curves with the fit parameters3
##Core genome
gamma = parameters_pan[0]
pan0 = parameters_pan[1]

genomes_curve = np.linspace(1, max(genomes))
pangenome_curve = pan0*(genomes_curve**gamma)
err_pan = np.sqrt(np.diag(parameters_cov_pan))


#Read data New genes
my_new_genes = pd.read_table("new_genes_micropan.csv",
                                   sep=",", 
                                   low_memory=False, header=None)
my_new_genes = my_new_genes.T
new_genes_melted = my_new_genes.melt()

genomes_new = new_genes_melted.iloc[:,0]
genomes_new = np.asarray(genomes_new)
genomes_new = genomes_new + 1
new_genes = new_genes_melted.iloc[:,1]

def new_genes_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_new, parameters_cov_new = curve_fit(new_genes_fit,
                                           genomes_new,
                                           new_genes, p0=[1, 1])

#Generate curves with the fit parameters3
##Core genome
alpha = parameters_new[0]
new0 = parameters_new[1]

genomes_curve_new = np.linspace(1, max(genomes_new))
new_genes_curve = new0*(genomes_curve_new**alpha)
err_core = np.sqrt(np.diag(parameters_cov_new))



#Plots
fig, ax = plt.subplots()
ax.plot(genomes, pangenome, '+', color = "bisque")
ax.plot(genomes_curve, pangenome_curve, 'darkorange', label = "Total")
ax.plot(genomes_new, new_genes, '+', color = "khaki")
ax.plot(genomes_curve, new_genes_curve, color = 'orange', label = "New")
ax.set_xlabel('Number of genomes', fontsize = 12)
ax.set_ylabel('Number of clusters', fontsize =12)

ax.text(1, 8000, r'$y = 2463.35 x ^ {0.27}$', fontsize = 12)
ax.text (1, 7000,
         r'$\gamma = {0:0.2f} \pm {1:0.4f}$'.format(gamma, err_pan[0]))

ax.text(60, 2000, r'$y = 569.97 x ^ {-0.69}$', fontsize = 12)
ax.text (60, 1000,
         r'$\alpha = {0:0.2f} \pm {1:0.3f}$'.format(-alpha, err_core[0]))
ax.grid()
ax.legend(loc="center right")

#plots 2
fig1, ((ax1, ax2)) = plt.subplots(ncols = 2, figsize = (7,3.5))
ax1.plot(genomes, pangenome, '+', color = "bisque")
ax1.plot(genomes_curve, pangenome_curve, 'darkorange', label = "Total")
ax2.plot(genomes_new, new_genes, '+', color = "khaki")
ax2.plot(genomes_curve, new_genes_curve, color = 'orange', label = "New")
ax1.set_xlabel('Number of genomes', fontsize= 12)
ax1.set_ylabel('Number of clusters', fontsize= 12)
ax1.set_xlim(1, 120)

ax2.set_xlabel('Number of genomes', fontsize= 12)
#ax2.set_ylabel('Number of clusters')
ax2.set_xlim(1, 120)

ax1.text(50, 5000, r'$y = 2463.35 x ^ {0.27}$', fontsize = 12)
ax1.text (50, 4200,
         r'$\gamma = {0:0.2f} \pm {1:0.4f}$'.format(gamma, err_pan[0]))

ax2.text(40, 500, r'$y = 569.97 x ^ {-0.69}$', fontsize = 12)
ax2.text (40, 400,
         r'$\alpha = {0:0.2f} \pm {1:0.3f}$'.format(-alpha, err_core[0]))
ax1.grid()
ax1.legend()
ax2.grid()
ax2.legend()

fig1.savefig("Micropan_new_total.svg", dpi=1200, format = 'svg')



