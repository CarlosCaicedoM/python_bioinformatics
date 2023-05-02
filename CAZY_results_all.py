# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 00:59:20 2021

@author: Carlos Caicedo-Montoya
"""

import pandas as pd
import re 
import numpy as np
import os

#my_directories
my_directories = [name for name in os.listdir(".") if os.path.isdir(name)]
CAZy_terms = []

for folder in my_directories:

    #read_data
    path = folder + "/overview.txt"
    my_data = pd.read_table(path, sep="\t")
    #Rename the column #ofTools by Num_of_Tools in order to avoid complications 
    #with the character #
    my_data.rename(columns={'#ofTools':'Num_of_Tools',}, 
                     inplace=True)
    #Filter the data to discard those annotations that not were predicted 
    #by more than 2 tools
    my_data_t = my_data[my_data['Num_of_Tools'] >= 2]
    
    #Eliminate the parenthesis and the data between them in the columns
    my_data_t['HMMER'] = my_data_t['HMMER'].map(lambda x: re.sub("[\(\[].*?[\)\]]", "", x))
    my_data_t['Hotpep'] = my_data_t['Hotpep'].map(lambda x: re.sub("[\(\[].*?[\)\]]", "", x))
    my_data_t['DIAMOND'] = my_data_t['DIAMOND'].map(lambda x: re.sub("[\(\[].*?[\)\]]", "", x))
    
    #Split the prediction of each different tools for each gene  
    #This is neccesary because some genes have more than one annotations in its 
    #sequence.
    #In addition, I eliminate all the sub-annotations 
    HMMER_expanded = my_data_t['HMMER'].str.split('+',expand=True)
    HMMER_expanded = HMMER_expanded.fillna(value='-')
    for i in range(HMMER_expanded.shape[1]):
        HMMER_expanded[i] = HMMER_expanded.iloc[:, i].map(lambda x: re.sub(r'_.*', "", x))
        
    
    Hotpep_expanded = my_data_t['Hotpep'].str.split('+',expand=True)
    Hotpep_expanded = Hotpep_expanded.fillna(value='-')
    for i in range(Hotpep_expanded.shape[1]):
        Hotpep_expanded[i] = Hotpep_expanded.iloc[:, i].map(lambda x: re.sub(r'_.*', "", x))
        
    
    DIAMOND_expanded = my_data_t['DIAMOND'].str.split('+',expand=True)
    DIAMOND_expanded = DIAMOND_expanded.fillna(value='-')
    for i in range(DIAMOND_expanded.shape[1]):
        DIAMOND_expanded[i] = DIAMOND_expanded.iloc[:, i].map(lambda x: re.sub(r'_.*', "", x))
        
    
    #Create a new dataframe with the predictions explite
    #Note that I discarded the SignalP results
    df = pd.concat([HMMER_expanded, Hotpep_expanded, DIAMOND_expanded], axis= 1)
    Gene_ID = my_data_t['Gene ID']
    df = df.set_index(Gene_ID)
    
    #Adjust the column names to the index between 0 and the shape of the dataframe
    columns = list(range(df.shape[1]))
    df.columns = columns
    
    df = df.T
    #Change the - values previously introduced by None values
    df = df.replace({'-': None})
    
    #Estimate the frequency of each annotation adding the predictions of different tools
    #This will allow to discard those domains predicted only by one software making the 
    #prediccions more accurate
    counts = []
    for label,  values in df.items():
        count = df[label].value_counts()
        counts.append(count)
    
    CAZy = pd.concat(counts, axis = 1)
    
    # This is an additional filter
    #all values equal to 1 in the frequency are discarded because they were 
    #predicted only by one software
    #Then the domains predicted by more than one tool or present more than two times
    #in a sequence are transformed to 1 (one) values to not carry out overestimation 
    #with the number of times a CAZy appears 
    max_value= int(max(list(CAZy.max())))
    my_results = CAZy.replace(1, 0)
    my_results = my_results.replace(list(range(2,max_value+1)), 1)
    my_results = my_results.replace(0, None)
    my_results[folder] = my_results.sum(axis=1)
    
    #Summarize the number of times the annotations appear in the genome
    my_results_def = my_results[my_results[folder] >= 1]
    my_results_def = my_results_def.loc[:, folder]
    CAZy_terms.append(my_results_def)
    


CAZy_matrix = pd.concat(CAZy_terms, axis = 1)

CAZy_matrix = CAZy_matrix.fillna(value=0)

#CAZy_matrix.to_csv("Cazy_Matrix.csv", index=True, header= True, sep=',')
 