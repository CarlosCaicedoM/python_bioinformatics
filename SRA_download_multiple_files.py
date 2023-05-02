#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 13:22:39 2022

@author: Carlos Caicedo-Montoya
"""

from IPython import get_ipython
get_ipython().magic('reset -sf')


import subprocess
import os
from pathlib import Path
import shutil


#Create file to download your data
my_experiment = "RNA-seq_Streptococcus"
Path_parent = Path(os.getcwd())
my_directory=str(Path_parent.parent) 
path = os.path.join(my_directory, "raw_data", my_experiment) 

try: 
    os.mkdir(path) 
except OSError as error: 
    print(error) 

my_file = "/SraAccList.txt"
my_data = my_directory+my_file
#If you want to remove all whitespace characters (newlines and spaces) from the end of each line
with open(my_data) as f:
    lines = [line.rstrip() for line in f]
Accessions  = lines[0:-1]

#Install SRA-tools
# conda install -c bioconda sra-tools 

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in Accessions:
    print ("Currently downloading: " + sra_id)
    fasterq = "fasterq-dump " + sra_id + " -O " + path
    print ("The command used was: " + fasterq)
    subprocess.call(fasterq, shell=True)


#Install FastQC
#conda install -c bioconda fastqc
#conda install -c bioconda/label/broken fastqc
#conda install -c bioconda/label/cf201901 fastqc
	
fastq_files = [file for file in os.listdir(path) if file.endswith(".fastq") ]
fastq_command =[]  
for i in fastq_files:
    fastq_command.append(os.path.join(path, i))
fastq_command.insert(0, "fastqc")
fastqc_outdir=os.path.join(my_directory, "results", "fastqc", my_experiment)

try:
    os.mkdir(os.path.join(my_directory, "results", "fastqc"))
except OSError as error:
    print(error)
try:
    os.mkdir(os.path.join(my_directory, "results", "fastqc", my_experiment))
except OSError as error:
    print(error)  
   
fastq_command.append("-o")
fastq_command.append(fastqc_outdir)
fastq_command.append("-t")
fastq_command.append("4")
#--noextract
subprocess.call(fastq_command)

#install multiQC
#conda install -c bioconda -c conda-forge multiqc

import multiqc
multiqc.run(fastqc_outdir)

try:
    os.mkdir(os.path.join(my_directory, "results", "multiqc"))
except OSError as error:
    print(error)
try:
    os.mkdir(os.path.join(my_directory, "results", "multiqc", my_experiment))
except OSError as error:
    print(error)  

multiQC_command = ["multiqc", fastqc_outdir]
multiQC_command.append("-n")
multiQC_command.append(my_experiment)
multiQC_command.append("-o")
multiQC_command.append(os.path.join(my_directory, "results", "multiqc", my_experiment))
multiQC_command.append("-f")
multiQC_command.append("--profile-runtime") 
subprocess.call(multiQC_command)

#Install cutadapt
#conda install -c bioconda cutadapt
#conda install -c bioconda/label/cf201901 cutadapt



fastq_files_path = []
for i in fastq_files:
    fastq_files_path.append(os.path.join(path, i))



try:
    os.mkdir(os.path.join(my_directory, "results", "cutadapt"))
except OSError as error:
    print(error)
try:
    os.mkdir(os.path.join(my_directory, "results", "cutadapt", my_experiment))
except OSError as error:
    print(error) 

adapter_5prime = "GTTCAGAGTTCTACAGTCCGACGATC"
adapter_3prime = "TGGAATTTCTCGGGTGCCAAGG"


#Cutadapt command for adapter trimming at both ends, 
#quality triming at both ends and filtering to only have reads with length >18nt


for i,j  in enumerate(fastq_files):
    cutadapt_command=["cutadapt", "-a", adapter_3prime, "-g",
                  adapter_5prime, "-q", "15,10", "-m", "18", 
                  "--cores=0", "-o",  
                  os.path.join(my_directory, "results", "cutadapt", my_experiment, "trimmed_"+j), fastq_files_path[i]]

    subprocess.call(cutadapt_command)
              
                 
                 
#Crear directorios raw_Data y results




#Simple quality filtering for FASTQ files
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

for sra_id in Accessions:
    fastq_file = sra_id+".fastq"
    count = 0
    #Count sequences
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq"):
        count += 1
    print("%i reads" % count)

    #Filter by quality 
    good_reads = (
        rec
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq")
        if min(rec.letter_annotations["phred_quality"]) >= 20
            )
    count = SeqIO.write(good_reads, "good_quality"+fastq_file, "fastq")
    print("Saved %i reads" % count)


#Trimming off primer sequences
count = 0
total_len = 0
with open(os.path.join(path, fastq_file)) as in_handle:
    for title, seq, qual in FastqGeneralIterator(in_handle):
        count += 1
        total_len += len(seq)
print("%i records with total sequence length %i" % (count, total_len))

#Sequences with the primers
primer_reads = (
    rec
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq")
    if rec.seq.startswith("GATGACGGTGT")
)
count = SeqIO.write(primer_reads, "with_primer.fastq", "fastq")
print("Saved %i reads" % count)


#Sequences with the primers but trimmed
trimmed_primer_reads = (rec[11:]
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq")
    if rec.seq.startswith("GATGACGGTGT"))
count = SeqIO.write(trimmed_primer_reads, "with_primer_trimmed.fastq", "fastq")
print("Saved %i reads" % count)



def trim_primer(record, primer):
    if record.seq.startswith(primer):
        return record[len(primer) :]
    else:
        return record


trimmed_reads = (trim_primer(record, "GATGACGGTGT")
    for record in SeqIO.parse(os.path.join(path, fastq_file), "fastq"))
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)


def trim_primers(records, primer):
    """Removes perfect primer sequences at start of reads.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    len_primer = len(primer)  # cache this for later
    for record in records:
        if record.seq.startswith(primer):
            yield record[len_primer:]
        else:
            yield record


original_reads = SeqIO.parse(os.path.join(path, fastq_file), "fastq")
trimmed_reads = trim_primers(original_reads, "GATGACGGTGT")
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)


#Trimming off adaptor sequences
#This time however, we will look for the sequence anywhere in the reads, 
#not just at the very beginning
def trim_adaptors(records, adaptor):
    """Trims perfect adaptor sequences.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    len_adaptor = len(adaptor)  # cache this for later
    for record in records:
        index = record.seq.find(adaptor)
        if index == -1:
            # adaptor not found, so won't trim
            yield record
        else:
            # trim off the adaptor
            yield record[index + len_adaptor :]
original_reads = SeqIO.parse(os.path.join(path, fastq_file), "fastq")
trimmed_reads = trim_adaptors(original_reads, "GATGACGGTGT")
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)

#let’s add a minimum length requirement as well:
def trim_adaptors(records, adaptor, min_len):
    """Trims perfect adaptor sequences, checks read length.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    len_adaptor = len(adaptor)  # cache this for later
    for record in records:
        len_record = len(record)  # cache this for later
        if len(record) < min_len:
            # Too short to keep
            continue
        index = record.seq.find(adaptor)
        if index == -1:
            # adaptor not found, so won't trim
            yield record
        elif len_record - index - len_adaptor >= min_len:
            # after trimming this will still be long enough
            yield record[index + len_adaptor :]

original_reads = SeqIO.parse("SRR020192.fastq", "fastq")
trimmed_reads = trim_adaptors(original_reads, "GATGACGGTGT", 100)
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)

