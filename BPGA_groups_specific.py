# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 23:47:27 2021

@author: Carlos Caicedo-Montoya
"""

#Import libraries

import numpy as np
import pandas as pd
from Bio import SeqIO


#Read data
my_data = pd.read_table("groups_specific_genes_3groups.txt", low_memory=False)
my_data.set_index("Status", inplace=True)

group1 = my_data.loc["PresentInGp_1"]

keep = ["GeneID", "Sequence"]
index_bool = np.isin(group1.columns, keep)
my_seqs = group1.iloc[:, index_bool]

#Save the table
my_seqs.to_csv('data.csv', sep='\t', header = False, index = False)

#Read with Biopython
records = SeqIO.parse("data.csv", "tab")
count = SeqIO.write(records, "plant_tissues.fasta", "fasta")
print("Converted %i records" % count)
