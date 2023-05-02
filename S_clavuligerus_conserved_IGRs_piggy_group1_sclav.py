#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 20:38:05 2021

@author: usuario
"""

#Import libraries
import pandas as pd 
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from os import remove
import subprocess
import numpy as np
import matplotlib
import seaborn as sns




#Read data
my_data_piggy_num = pd.read_table('IGR_presence_absence.Rtab',
                     sep='\t',
                     low_memory=False)
# Change dataframe index
my_data_piggy_num.set_index('Gene', inplace=True)

#Establish the 
addition = my_data_piggy_num.sum(axis=1) 
my_data_piggy_num['addition'] = addition


my_group = {"S._lunaelactis_MM109":"GCF_003054555.1_ASM305455v1_genomic.gbff", 
            "S._pristinaespiralis_HCCB_10218":"GCF_001278075.1_ASM127807v1_genomic.gbff",
            "S._clavuligerus_ATCC_27064":"GCF_005519465.1_ASM551946v1_genomic.gbff",
            "S._clavuligerus_F1D-5":"GCF_003454755.1_ASM345475v1_genomic.gbff",
            "S._clavuligerus_F613-1":"GCF_001693675.1_ASM169367v1_genomic.gbff"}

my_group2 = {y:x for x,y in my_group.items()}
my_assemblies = list(my_group.values())




#Rename the columns  
#This  is important to select the cluster only present in the genome
#of S.clavuligerus ATCC 27064

my_data_piggy_num.rename(columns = my_group2, inplace = True)

#Select the IGR clusters present only in Sclavuligerus ATCC27064
my_IGRs = my_data_piggy_num[my_data_piggy_num["S._clavuligerus_ATCC_27064"] == 1]

#Select the IGRs conserved at least in 4 genomes
my_conserved_IGRs = my_IGRs[my_IGRs["addition"] >=4]



#Read data with the information of sequences names
my_data_piggy = pd.read_table('IGR_presence_absence.csv',
                     sep=',',
                     low_memory=False)

# Set index (group name)
my_data_piggy.set_index('Gene', inplace=True)

# Drop the other info columns
my_data_piggy.drop(list(my_data_piggy.columns[:13]), axis=1, inplace=True)
#my_matrix.drop('addition', axis=1, inplace=True)

#Muscle must be installed in your current directory
muscle_exe = "./muscle"

Lenghts =[]

for i in range(my_conserved_IGRs.shape[0]):
    IGRs_total = my_data_piggy.iloc[i]
    IGRs = IGRs_total.tolist()
    file = my_data_piggy.index[i]

    with open("headers_"+file + ".txt", 'w') as output:
        for row in IGRs:
            output.write(str(row) + '\n')
    

    #Extract IGR clusters

    input_file = "IGR_sequences.fasta"
    id_file = "headers_"+ file + ".txt"
    output_file = file + ".fasta"
   
    with open(id_file) as id_handle:
        wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
    print("Found %i unique identifiers in %s" % (len(wanted), id_file))

    #+_+_ is used as a delimiter between fields in piggy
    records = (r for r in SeqIO.parse(input_file, "fasta") if r.description in wanted)
    
    count = SeqIO.write(records, output_file, "fasta")
    print("Saved %i records from %s to %s" % (count, input_file, output_file))
    if count < len(wanted):
        print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))

    x = list(len(r) for r in SeqIO.parse(output_file, "fasta"))
    Lenghts.append(x[0])
    if x[0] >= 50:
        muscle_cline = MuscleCommandline(muscle_exe, input=output_file, 
                                     out="muscle_" + output_file + ".aln",
                                     clwstrict=True)
        stdout, stderr = muscle_cline()

    
    remove("headers_"+file + ".txt")
    remove(output_file)

Lenghts =np.asarray(Lenghts)
sns.displot(Lenghts)
min(Lenghts)
max(Lenghts)

subprocess.call('./run_RNAz_local.sh')







