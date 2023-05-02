#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 12:52:17 2022

@author: Carlos Caicedo-Montoya
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

#%% Extract rRNAs from gen bank

from IPython import get_ipython
get_ipython().magic('reset -sf')


import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os
import subprocess



handle = "/media/usuario/lab_bioinformatica/Genome-wide_identification_sRNAS_SClavuligerus/raw_data/Genbank_group1/GCF_003054555.1_ASM305455v1_genomic.gbff"

rRNAs = []
    
for seq_record in SeqIO.parse(handle, "genbank"):
    for feature in seq_record.features:
        if feature.type == 'rRNA':
            my_start = feature.location.start.position
            my_end = feature.location.end.position
            gene_ID = str(feature.qualifiers['locus_tag'][0] + " " + feature.qualifiers['product'][0])                       
            rRNA_seq = seq_record.seq[my_start:my_end]
            rRNAs.append(SeqRecord(rRNA_seq,
                                        id="%s" % (gene_ID), 
                                        description=""))
            
            
            
base_name=seq_record.dbxrefs[2]
base_name_list=base_name.split(":")
base_name_formated = base_name_list[1]
outpath = os.path.splitext(os.path.basename(base_name_formated))[0] + "_RNAs.fasta"
SeqIO.write(rRNAs, open(outpath,"w"), "fasta")            
 
            
import subprocess 
subprocess.call(['grep', '-c', '>', 'GCF_003054555_RNAs.fasta'])
           