# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 22:21:40 2020

@author:  Carlos Caicedo-Montoya
"""

#Import libraries
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Read data
my_data = pd.read_table('Core-Pan-Genome-Streptomyces.txt',
                     sep='\t',
                     low_memory=False)

number_of_genomes = my_data.iloc[:, 0 ]
number_of_genomes = np.asarray(number_of_genomes)

number_of_genomes2 = my_data.iloc[:, 2 ]
number_of_genomes2 = np.asarray(number_of_genomes2)

pan_genome =  my_data.iloc[:, 1]
pan_genome = np.asarray(pan_genome)
core_genome = my_data.iloc[:, 3]
core_genome = np.asarray(core_genome)

#Fit power law curve
#PANGENOME

def pangenome_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_pan, parameters_cov_pan = curve_fit(pangenome_fit,
                                           number_of_genomes,
                                           pan_genome, p0=[1, 1])

def coregenome_fit(x, alpha, G0):
    return G0*(x**alpha)
parameters_core, parameters_cov_core = curve_fit(coregenome_fit,
                                           number_of_genomes,
                                           core_genome, p0=[-0.4, 3400])

#Generate curves with the fit parameters3
##Pangenome
gamma = parameters_pan[0]
pan0 = parameters_pan[1]
genomes = np.linspace(1, max(number_of_genomes))
pangenome_curve = pan0*(genomes**gamma)
err_pan = np.sqrt(np.diag(parameters_cov_pan))

##Core genome
alpha = parameters_core[0]
core0 = parameters_core[1]

genomes = np.linspace(1, max(number_of_genomes))
coregenome_curve = core0*(genomes**alpha)
err_core = np.sqrt(np.diag(parameters_cov_core))

#Plots
fig, ax = plt.subplots(figsize = (8, 6))
ax.plot(number_of_genomes, pan_genome, 'dy')
ax.plot(number_of_genomes, core_genome, 'dc')
ax.plot(genomes, pangenome_curve, 'darkred', label = 'Pangenome power law fit')
ax.plot(genomes, coregenome_curve, 'darkgreen', label = "Core genome power law fit")
ax.set_xlabel('Number of genomes', fontsize = 12)
ax.set_ylabel('Number of gene families', fontsize =12)
ax.text(20, 5500, r'$y = 5299.36 x ^ {-0.62}$', fontsize = 12)
ax.text(20, 30000, r'$y = 8923.88 x ^ {0.52}$', fontsize = 12)
ax.text (22, 3500,
         r'$\alpha = {0:0.2f} \pm {1:0.3f}$'.format(alpha, err_core[0]))
ax.text (22, 28000,
         r'$\gamma = {0:0.2f} \pm {1:0.3f}$'.format(gamma, err_pan[0]))


ax.grid()
ax.legend()

#Save
fig.savefig("core-pan-streptomyces", dpi =1200, format = 'svg')
