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
my_total_genes = pd.read_table("number_of_genes_in_pan_genome.txt",
                                   sep="\t", 
                                   low_memory=False, header=None)
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
my_new_genes = pd.read_table("number_of_new_genes.txt",
                                   sep="\t", 
                                   low_memory=False, header=None)
new_genes_melted = my_new_genes.melt()

genomes = new_genes_melted.iloc[:,0]
genomes = np.asarray(genomes)
genomes = genomes + 1
new_genes = new_genes_melted.iloc[:,1]

def new_genes_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_new, parameters_cov_new = curve_fit(new_genes_fit,
                                           genomes,
                                           new_genes, p0=[1, 1])

#Generate curves with the fit parameters3
##Core genome
alpha = parameters_new[0]
new0 = parameters_new[1]

genomes_curve = np.linspace(1, max(genomes))
new_genes_curve = new0*(genomes_curve**alpha)
err_core = np.sqrt(np.diag(parameters_cov_new))



#Plots
fig, ax = plt.subplots()
ax.plot(genomes, pangenome, '+', color = "lightsteelblue")
ax.plot(genomes_curve, pangenome_curve, 'royalblue', label = "Total genes")
ax.plot(genomes, new_genes, '+', color = "lightsteelblue")
ax.plot(genomes_curve, new_genes_curve, color = 'deepskyblue', label = "New genes")
ax.set_xlabel('Number of genomes', fontsize = 12)
ax.set_ylabel('Number of genes', fontsize =12)

ax.text(1, 120000, r'$y = 8212.43 x ^ {0.60}$', fontsize = 12)
ax.text (1, 110000,
         r'$\gamma = {0:0.2f} \pm {1:0.3f}$'.format(gamma, err_pan[0]))

ax.text(60, 20000, r'$y = 5993.29 x ^ {-0.45}$', fontsize = 12)
ax.text (60, 6000,
         r'$\alpha = {0:0.2f} \pm {1:0.3f}$'.format(-alpha, err_core[0]))
ax.grid()
ax.legend(loc="center right")

#plots 2
fig1, ((ax1, ax2)) = plt.subplots(ncols = 2, figsize = (7,3.5))
ax1.plot(genomes, pangenome, '+', color = "lightsteelblue")
ax1.plot(genomes_curve, pangenome_curve, 'royalblue', label = "Total")
ax2.plot(genomes, new_genes, '+', color = "lightsteelblue")
ax2.plot(genomes_curve, new_genes_curve, color = 'deepskyblue', label = "New")
ax1.set_xlabel('Number of genomes', fontsize= 12)
ax1.set_ylabel('Number of clusters', fontsize= 12)
ax1.set_xlim(1, 120)

ax2.set_xlabel('Number of genomes', fontsize= 12)
#ax2.set_ylabel('Number of clusters')

ax1.text(50, 40000, r'$y = 8212.43 x ^ {0.60}$', fontsize = 12)
ax1.text (50, 30000,
         r'$\gamma = {0:0.2f} \pm {1:0.3f}$'.format(gamma, err_pan[0]))

ax2.text(40, 6500, r'$y = 5993.29 x ^ {-0.45}$', fontsize = 12)
ax2.text (40, 5800,
         r'$\alpha = {0:0.2f} \pm {1:0.3f}$'.format(-alpha, err_core[0]))
ax1.grid()
ax1.legend()
ax2.grid()
ax2.legend()
ax2.set_xlim(1, 120)

fig1.savefig("Roary_new_total.svg", dpi=1200, format = 'svg')




