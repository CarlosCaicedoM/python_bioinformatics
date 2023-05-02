#!/usr/bin/env python3
"""
Created on Wed Aug 11 17:41:07 2021

@author: Carlos Caicedo-Montoya

@author: Carlos Caicedo-Montoya
"""

#Import libraries
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns



#Read data
my_data_piggy = my_data_piggy = pd.read_csv('IGR_presence_absence.Rtab_all_pangenome.csv')


my_group = {"GID121":"S._clavuligerus_ATCC_27064", 
            "GID101":"S._clavuligerus_F1D-5", 
            "GID46":"S._clavuligerus_F613-1",
            "GID83":"S._lunaelactis_MM109",
            "GID34":"S._pristinaespiralis_HCCB_10218"}


my_group2 = {y:x for x,y in my_group.items()}



my_assemblies = list(my_group2.values())


# Change dataframe index
my_data_piggy.set_index("Gene", inplace=True)

#Select the desired genomes
my_group_piggy = my_data_piggy[my_data_piggy.columns.intersection(my_assemblies)]

my_group_piggy2 = my_data_piggy[[c for c in my_data_piggy.columns if c in my_assemblies]]

my_group_piggy.equals(my_group_piggy2)

my_group_piggy.rename(columns = my_group, inplace = True)


my_group_piggy["addition"]=my_group_piggy.sum(axis = 1)

my_IGRs = my_group_piggy[my_group_piggy["addition"] >= 1]

my_IGRs = my_IGRs[my_IGRs["S._clavuligerus_ATCC_27064"] == 1]

my_conserved_IGRs = my_IGRs[my_IGRs["addition"] >= 2]

my_conserved_IGRs2 = my_IGRs[my_IGRs["addition"] >= 3]

my_conserved_IGRs3 = my_IGRs[my_IGRs["addition"] >= 4]

my_conserved_IGRs4 = my_IGRs[my_IGRs["addition"] >=5]









