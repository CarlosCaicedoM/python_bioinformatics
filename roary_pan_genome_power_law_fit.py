# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 19:22:20 2021

@author: Carlos Caicedo-Montoya
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

#Read data

my_conserved_genes = pd.read_table("number_of_genes_in_pan_genome.txt",
                                   sep="\t", 
                                   low_memory=False, header=None)
df_melted = my_conserved_genes.melt()


genomes = df_melted.iloc[:,0]
genomes = np.asarray(genomes)
genomes = genomes + 1
genes = df_melted.iloc[:,1]


#Fit power law curve
#PANGENOME

def pangenome_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_pan, parameters_cov_pan = curve_fit(pangenome_fit,
                                           genomes,
                                           genes, p0=[1, 1])

#Generate curves with the fit parameters3
##Core genome
alpha = parameters_pan[0]
pan0 = parameters_pan[1]

genomes_curve = np.linspace(1, max(genomes))
pangenome_curve = pan0*(genomes_curve**alpha)
err_core = np.sqrt(np.diag(parameters_cov_pan))


#Plots
fig, ax = plt.subplots()
ax.plot(genomes, genes, '+g')
ax.plot(genomes_curve, pangenome_curve, 'maroon', label = "Total genes")
ax.set_xlabel('Number of genomes', fontsize = 12)
ax.set_ylabel('Number of genes', fontsize =12)
ax.text(60, 40000, r'$y = 8212.43 x ^ {0.60}$', fontsize = 12)
ax.text (60, 30000,
         r'$\gamma = {0:0.2f} \pm {1:0.3f}$'.format(alpha, err_core[0]))
ax.grid()
ax.legend()


