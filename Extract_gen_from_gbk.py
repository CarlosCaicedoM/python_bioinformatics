# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 03:48:00 2022

@author: Carlos Caicedo-Montoya
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


handle = "D:\GCF_005519465.1_ASM551946v1_genomic.gbk"
record_iterator = SeqIO.parse(handle, "genbank")
first_record = next(record_iterator)

region_AC = first_record.seq[4800000:4900000]




region_CA = SeqRecord(region_AC, 
                       id="region_CA",
                       description="")

SeqIO.write(region_CA, "region_CA_old.fasta", "fasta")            
