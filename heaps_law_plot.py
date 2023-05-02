# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 00:19:52 2021

@author: Carlos Caicedo-Montoya
"""
import numpy as np
import matplotlib.pyplot as plt

#Roary
k_r = 3967.51
alpha_r = 0.3284
N = np.arange(1,122)

n = k_r*(N**(-alpha_r))

#Micropan 

k_m = 882.238
alpha_m = 0.8153
N = np.arange(1,122)

n_m= k_m*(N**(-alpha_m))


#fig
fig, ax = plt.subplots()
ax.plot(N, n, label = "Roary")
ax.plot(N, n_m, color = "darkorange", label = "Micropan")
ax.set_xlabel("Number of genomes", fontsize=10)
ax.set_ylabel("Number of new genes", fontsize=10)
ax.grid()
ax.legend()