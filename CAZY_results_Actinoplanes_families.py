# -*- coding: utf-8 -*-
"""
Created on Fri May 21 09:46:37 2021

@author: Carlos Caicedo-Montoya
"""
# %% 
import pandas as pd
import re 
import numpy as np
import os
from scipy.stats import ranksums

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
    #In addition, I will eliminate all the sub-annotations 
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
    #my_results = my_results.replace(0, None)
    my_results = my_results.replace({0:None})
    my_results[folder] = my_results.sum(axis=1)
    
    #Summarize the number of times the annotations appear in the genome
    my_results_def = my_results[my_results[folder] >= 1]
    my_results_def = my_results_def.loc[:, folder]
    CAZy_terms.append(my_results_def)
    


CAZy_matrix = pd.concat(CAZy_terms, axis = 1)

CAZy_matrix = CAZy_matrix.fillna(value=0)

#CAZy_matrix.to_csv("Cazy_Matrix_all.csv", index=True, header= True, sep=',')


#Define which genomes are from endophytes strains 
endophytes = ["Actinoplanes.sp.TFC3_out_CGC"]
others = []
for i in my_directories:
    if i not in endophytes:
        others.append(i)

#Families

terms = CAZy_matrix.index
families = ["AA", "CBM", "CE", "GH", "GT", "PL"]
wilk_test = []
##AA
for family in families:
    #Extract the terms belonging to an specified family using 
    #list comprehension
    family_terms = [x for x in terms if family in x]
    CAZy_family = CAZy_matrix.loc[family_terms, :]
    #Add the values for each family in each genome
    total_family = CAZy_family.sum()
    total_family = total_family.to_dict()
    total_family_endophytes = []
    total_family_others = []
    #Assign each count to a list of endophytes or no endophytes 
    for i in total_family.keys():
        if i in endophytes:
            total_family_endophytes.append(total_family[i])
        else:
            total_family_others.append(total_family[i])
    #Calculate the 
    wilk_test.append(ranksums(total_family_endophytes, total_family_others))


for i, j in enumerate(wilk_test):
    print (families[i], ",  p_value =", j[1])
    

# %%Plots
count_families = []
for family in families:
    #Extract the terms belonging to an specified family using 
    #list comprehension
    family_terms = [x for x in terms if family in x]
    CAZy_family = CAZy_matrix.loc[family_terms, :]
    #Add the values for each family in each genome
    total_family = CAZy_family.sum()
    count_families.append(total_family)
    

CF = pd.concat(count_families, axis=1)  
CF = CF.set_axis(families, axis="columns")

strains = ["TFC3", "SE50/110", "431", "N902-109", "7358",
           "SE50", "ATCC 31121", "OR16", "SE50/110 ACP50"]


CF_copy = CF.copy()
CF_copy = CF_copy.set_axis(strains, axis = "index")
CF_copy = CF_copy.astype(int)


import seaborn as sns
import matplotlib 
import matplotlib.pyplot as plt
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)
fig, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(CF_copy, cmap="PiYG", annot=True, fmt = "d", ax = ax,
            linewidths=.8) 

#fig.savefig("summary_CAZy_Actinoplanes.pdf", format = "pdf", dpi = 1200)


# %% PCWDE, plant cell wall degrading enzymes

chitin = ["GH18", "GH19", "GH48", "GH23"]
cellulose  = ["GH5", "GH6", "GH7",  "GH8", "GH9", "GH45", "AA10"]
hemicellulose = ["GH10", "GH11", "GH12", "GH16", "GH29", "GH30", "GH39", "GH62", "GH67", 
         "GH74", "GH95", "GH115", "GH120", "GH141", "GH142", "CE1", "CE2", 
         "CE3", "CE5", "CE7", "CE15"]
galactomannan = ["GH26", "GH27", "GH36", "GH76", "GH92", "GH113", "GH130", "GH134"]
starch = ["GH13", "GH15", "GH31", "GH77", "GH133"]
pectin = ["GH1", "GH2", "GH28", "GH42", "GH43", "GH53", "GH78", "GH88", "GH93",
          "GH105", "PL1", "PL3", "PL4", "PL9", "PL11", "PL22", "CE8", "CE12"]
inulin = ["GH32"]
peptidoglycan = ["GH3", "GH25"]



PP = {"chitin":["GH18", "GH19", "GH48", "GH23"], 
      "cellulose":   ["GH5", "GH6",  "GH7",  "GH8", "GH9", "GH45", "AA10"], 
      "hemicellulose" : ["GH10", "GH11", "GH12", "GH16", "GH29", "GH30", "GH39", 
                         "GH62", "GH67", "GH74", "GH95", "GH115", "GH120", 
                         "GH141", "GH142", "CE1", "CE2", "CE3", "CE5", "CE7",
                         "CE15"], 
      "galactomannan" :["GH26", "GH27", "GH36", "GH76", "GH92", "GH113", 
                        "GH130", "GH134"], 
      "starch" : ["GH13", "GH15", "GH31", "GH77", "GH133"], 
      "pectin" : ["GH1", "GH2", "GH28", "GH42", "GH43", "GH53", "GH78",
                  "GH88", "GH93", "GH105", "PL1", "PL3", "PL4", "PL9", "PL11",
                  "PL22", "CE8", "CE12"], 
      "inulin" : ["GH32"], 
      "peptidoglycan" : ["GH3", "GH25"]}

PP_summary ={}
##AA
for key, value in PP.items():
    #Extract the terms belonging to an specified family using 
    #list comprehension
    test=[]
    for family in value:
        PP_terms = [x for x in terms if family == x]
        CAZy_PP = CAZy_matrix.loc[PP_terms, :]
        test.append(CAZy_PP)
        result = pd.concat(test)
        result= result.sum()
    PP_summary[key] = result
        

PP_final = pd.concat(PP_summary)
my_matrix = PP_final.unstack()

strains = ["TFC3", "SE50/110", "431", "N902-109", "7358",
           "SE50", "ATCC 31121", "OR16", "SE50/110 ACP50"]


PP_copy = my_matrix.copy()
PP_copy = PP_copy.T
PP_copy = PP_copy.set_axis(strains, axis = "index")
PP_copy = PP_copy.astype(int)


import seaborn as sns
import matplotlib 
import matplotlib.pyplot as plt
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)
fig, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(PP_copy.T, cmap="viridis", annot=True, fmt = "d", ax = ax,
            linewidths=.8) 