# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 22:55:18 2021

@author: Carlos Caicedo-Montoya
"""
import pandas as pd
import numpy as np


#read data

ANI_matrix = pd.read_table("ANI_matrix.csv",
                                   sep=",", 
                                   low_memory=False)

ANI_matrix.set_index("reference", inplace=True)
ANI_matrix = ANI_matrix.where(np.tril(np.ones(ANI_matrix.shape), k=-1).astype(np.bool))
ANI_matrix = ANI_matrix/100

ANI_matrix.to_csv("lower_ANI_Matrix", sep=',', header =True, index=True)