# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 23:27:29 2021

@author: User
"""

from Bio import SeqIO

input_file = "representative_sequences.fa"
id_file = "cloud_headers.txt"
output_file = "cloud_proteins"

with open(id_file) as id_handle:
    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
print("Found %i unique identifiers in %s" % (len(wanted), id_file))

records = (r for r in SeqIO.parse(input_file, "fasta") if r.description.split()[0] in wanted)
count = SeqIO.write(records, output_file, "fasta")
print("Saved %i records from %s to %s" % (count, input_file, output_file))
if count < len(wanted):
    print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))
    
    
