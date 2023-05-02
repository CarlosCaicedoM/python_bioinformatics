# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 22:33:41 2020

@author: Carlos Caicedo-Montoya
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Convert FASTA to tab (or other formats)

records = SeqIO.parse("core_gene_alignment.aln.fasta", "fasta")
count = SeqIO.write(records, "core_gene_alignment_names.aln.tab", "tab")
print("Converted %i records" % count)


#Change the header of a fasta file

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

my_data = pd.read_table('strain_names.txt',
                     sep='\t',
                     low_memory=False)
strains = my_data['strains'].to_list()

my_sequences = []
for i, seq_record in enumerate(SeqIO.parse("core_gene_alignment.aln.fasta", "fasta")):
    rec = SeqRecord(seq_record.seq, id = strains[i], description="")
    my_sequences.append(rec)

count = SeqIO.write(my_sequences, "core_gene_alignment_names.aln.fasta", "fasta")
print("Written %i records" % count)



#Number of nucleotides in a fasta file
from Bio import SeqIO
print(sum(len(r) for r in SeqIO.parse("pan_genome_reference.fa", "fasta")))

#Number of sequences in a fasta file
from Bio import SeqIO
my_dict = SeqIO.index("pan_genome_reference.fa", "fasta")
print(len(my_dict))




