# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 15:49:49 2021
@author: Carlos Caicedo-Montoya
"""

#Import libraries
import pandas as pd

#Read data
my_data = pd.read_table("core_micropan_sep.txt", sep = '\t', low_memory =False)
my_data.set_index("Protein_Accession", inplace=True)

Protein_Accession = my_data.index
Protein_Accession = list(Protein_Accession)

sequences=[]
for i in Protein_Accession:
    if i not in sequences:
        sequences.append(i)

my_terms={}
for sequence in sequences:
    df1 = my_data.groupby(my_data.index).get_group(sequence)
    df2 = df1.stack(level=0)
    GO_terms = df2.tolist()
    true_terms = []
    for i in GO_terms:
        if i not in true_terms:
            true_terms.append(i)
    my_terms[sequence] = true_terms


my_terms_df=pd.DataFrame.from_dict(my_terms, orient='index')

my_terms_df.to_csv("GO_terms_core_micropan.csv", index=True, sep=',')