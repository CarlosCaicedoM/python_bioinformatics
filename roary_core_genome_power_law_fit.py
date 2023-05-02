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

my_conserved_genes = pd.read_table("number_of_conserved_genes.Rtab.txt",
                                   sep="\t", 
                                   low_memory=False, header=None)
df_melted = my_conserved_genes.melt()


genomes = df_melted.iloc[:,0]
genomes = np.asarray(genomes)
genomes = genomes + 1
genes = df_melted.iloc[:,1]


#Fit power law curve
#PANGENOME

def coregenome_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_core, parameters_cov_core = curve_fit(coregenome_fit,
                                           genomes,
                                           genes, p0=[1, 1])

#Generate curves with the fit parameters3
##Core genome
alpha = parameters_core[0]
core0 = parameters_core[1]

genomes_curve = np.linspace(1, max(genomes))
coregenome_curve = core0*(genomes_curve**alpha)
err_core = np.sqrt(np.diag(parameters_cov_core))


#Plots
fig, ax = plt.subplots()
ax.plot(genomes, genes, 'dc')
ax.plot(genomes_curve, coregenome_curve, 'darkgreen', label = "Core genome power law fit")
ax.set_xlabel('Number of genomes', fontsize = 12)
ax.set_ylabel('Number of genes', fontsize =12)
ax.text(40, 4000, r'$y = 4612.06 x ^ {-0.39}$', fontsize = 12)
ax.text (40, 3500,
         r'$\alpha = {0:0.2f} \pm {1:0.3f}$'.format(alpha, err_core[0]))
ax.grid()
ax.legend()
