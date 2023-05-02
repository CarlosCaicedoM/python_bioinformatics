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

my_conserved_genes = pd.read_table("number_of_new_genes.txt",
                                   sep="\t", 
                                   low_memory=False, header=None)
df_melted = my_conserved_genes.melt()


genomes = df_melted.iloc[:,0]
genomes = np.asarray(genomes)
genomes = genomes + 1
genes = df_melted.iloc[:,1]


#Fit power law curve
#PANGENOME

def new_genes_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_new, parameters_cov_new = curve_fit(new_genes_fit,
                                           genomes,
                                           genes, p0=[1, 1])

#Generate curves with the fit parameters3
##Core genome
alpha = parameters_new[0]
new0 = parameters_new[1]

genomes_curve = np.linspace(1, max(genomes))
new_genes_curve = new0*(genomes_curve**alpha)
err_core = np.sqrt(np.diag(parameters_cov_new))


#Plots
fig, ax = plt.subplots()
ax.plot(genomes, genes, '+y')
ax.plot(genomes_curve, new_genes_curve, 'darkmagenta', label = "New genes")
ax.set_xlabel('Number of genomes', fontsize = 12)
ax.set_ylabel('Number of genes', fontsize =12)
ax.text(60, 7000, r'$y = 5993.29 x ^ {-0.45}$', fontsize = 12)
ax.text (60, 6000,
         r'$\alpha = {0:0.2f} \pm {1:0.3f}$'.format(alpha, err_core[0]))
ax.grid()
ax.legend()


