# -*- coding: utf-8 -*-

"""
Created on Thu Jan 21 00:02:48 2021
@author: Carlos Caicedo-Montoya
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

#Read data

my_conserved_genes = pd.read_table("number_of_unique_genes.txt",
                                   sep="\t", 
                                   low_memory=False, header=None)
df_melted = my_conserved_genes.melt()


genomes = df_melted.iloc[:,0]
genomes = np.asarray(genomes)
genomes = genomes + 1
genes = df_melted.iloc[:,1]


#Fit power law curve
#PANGENOME

def unique_genes_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_unique, parameters_cov_unique = curve_fit(unique_genes_fit,
                                           genomes,
                                           genes, p0=[1, 1])

#Generate curves with the fit parameters3
##Core genome
alpha = parameters_unique[0]
unique0 = parameters_unique[1]

genomes_curve = np.linspace(1, max(genomes))
unique_genes_curve = unique0*(genomes_curve**alpha)
err_core = np.sqrt(np.diag(parameters_cov_unique))


#Plots
fig, ax = plt.subplots()
ax.plot(genomes, genes, '+b', linewidth = 1.5)
ax.plot(genomes_curve, unique_genes_curve, 'orange', 
        label = "Unique genes", linewidth = 2.0)
ax.set_xlabel('Number of genomes', fontsize = 12)
ax.set_ylabel('Number of genes', fontsize =12)
ax.text(60, 30000, r'$y = 6937.95 x ^ {0.51}$', fontsize = 12)
ax.text (60, 20000,
         r'$\gamma = {0:0.2f} \pm {1:0.3f}$'.format(alpha, err_core[0]))
ax.grid()
ax.legend()


