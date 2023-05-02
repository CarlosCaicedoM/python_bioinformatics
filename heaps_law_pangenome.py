# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 13:46:00 2020

@author: Carlos Caicedo-Montoya

A simple explanation of the power law for pangenome size characterization
"""

import numpy as np
import matplotlib.pyplot as plt

#Number of gene families
gamma = np.linspace(0.01, 0.3, 6)
alpha = 1-gamma
N = list(range(1,101))
N=np.asarray(N)
k=7811
n_array =np.zeros(len(gamma)*len(N)).reshape(len(gamma), len(N))

for i, j in enumerate(gamma):
    n_array[i] = k*N**(j)

#Figures
fig, ax = plt.subplots()
for i in range(len(n_array)):
    ax.plot(N, n_array[i], label=r'$ \gamma = {0:0.2f}$'.format(gamma[i]))
ax.set_xlabel('Number of genomes')
ax.set_ylabel('Numbero of gene families')
ax.legend(loc='upper left')
ax.grid()


#Number of new gene families
alpha = np.linspace(0.01, 2, 6)
N = list(range(1,101))
N=np.asarray(N)
k=4200
n_array =np.zeros(len(alpha)*len(N)).reshape(len(alpha), len(N))

for i, j in enumerate(alpha):
    n_array[i] = k*N**(-j)

#Figures
fig, ax = plt.subplots()
for i in range(len(n_array)):
    ax.plot(N, n_array[i], label=r'$ \alpha = {0:0.2f}$'.format(alpha[i]))
ax.set_xlabel('Number of genomes')
ax.set_ylabel('Numbero of new gene families')
ax.legend(loc='center')
ax.grid()








