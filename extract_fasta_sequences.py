# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 23:27:29 2021

@author: Carlos Caicedo-Montoya

Adapted from Byopython tutorial available at:
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec379
"""

from Bio import SeqIO

input_file = "pangenome_reference.fa"
id_file = "headers.txt"
output_file = "my_sequences.fa"

with open(id_file) as id_handle:
    wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
print("Found %i unique identifiers in %s" % (len(wanted), id_file))

records = (r for r in SeqIO.parse(input_file, "fasta") if r.description.split()[0] in wanted)
count = SeqIO.write(records, output_file, "fasta")
print("Saved %i records from %s to %s" % (count, input_file, output_file))
if count < len(wanted):
    print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))
    
    
