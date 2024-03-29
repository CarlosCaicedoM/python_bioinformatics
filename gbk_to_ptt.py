# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 07:46:42 2022

@author: Carlos Caicedo-Montoya
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')


"""
Title: GenBank file to Protein Table file (.ptt) and RNA table file (.rnt) parser
Input files: [./sequence.gb, ./*.gb (batch)]

This script creates a Protein Table file (.ptt) and a RNA table file (.rnt) from the given GenBank file
Multiple files can be given (or using *) for batch processing

@author: Henry Wiersma UMCG, Groningen
@date: 7-11-2017
@version: 1.0
"""

from Bio import SeqIO
import os
import sys
from random import randint


#rna feature types in the GenBank file (https://www.ncbi.nlm.nih.gov/books/NBK293913/)
rntFeatureTypes = ["rna", "mRNA", "tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA"]

#Template of the header of the tables
pttHeader="{description} - 0..{length}\n{numRows} proteins\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n"
rntHeader="{description} - 0..{length}\n{numRows} RNAs\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n"

#template of the rows of the table
#row="{location}\t{strand}\t{length}\t{pid}\t{gene}\t{synonym}\t{code}\t{cog}\t{product}\r\n"
row="{location}\t{strand}\t{length}\t{pid}\t{gene}\t{synonym}\t{code}\t{cog}\t{product}"


#createFile(pttOutputFilePath, pttHeader, description, length, pttRows)
#filePath = pttOutputFilePath
#headerTemplate=pttHeader
#rows = pttRows


def createFile(filePath, headerTemplate, description, length, rows):
    outputFile = open(filePath, 'w')
    header = headerTemplate.format(description=description,
                                   length=length,
                                   numRows=len(rows))
    outputFile.write(header)
    for row in rows:
        outputFile.write("%s\n" % row)
    outputFile.close()
    


def processFile(filePath, outputPath):
    file = SeqIO.parse(filePath, "genbank")
    fileName = os.path.basename(filePath)
    fileNameBase = fileName.split(".")[0]
    gbItem = 0
    print("Process file", fileName)
    for seqRec in file:
        pttRows = []
        rntRows = []

        # get description of the gb file
        description = "Description not available"
        if seqRec.description and seqRec.description != "":
            description = seqRec.description
        elif seqRec.id and seqRec.id != "":
            description = seqRec.id

        # get the sequence length
        length = 0
        if seqRec.seq:
            length = len(seqRec.seq)
        a =  randint(100000000, 500000000)     # randint is inclusive at both ends
        pid = a
        for i, feature in enumerate(seqRec.features):
            # location
            paramLocation = "{}..{}".format((feature.location.start + 1), feature.location.end)

            # strand
            paramStrand = "+" if feature.location.strand == 1 else "-"

            # gene name
            paramGene = "-"
            if "gene" in feature.qualifiers and len(feature.qualifiers["gene"]) > 0 and feature.qualifiers["gene"][
                0] != "":
                paramGene = (feature.qualifiers["gene"][0])

            # locus tag (Synonym)
            paramSynonym = "-"
            if "locus_tag" in feature.qualifiers and len(feature.qualifiers["locus_tag"]) > 0 and \
                            feature.qualifiers["locus_tag"][0] != "":
                paramSynonym = (feature.qualifiers["locus_tag"][0])

            # product descriptopn
            paramProduct = "-"
            if "product" in feature.qualifiers and len(feature.qualifiers["product"]) > 0 and \
                            feature.qualifiers["product"][
                                0] != "":
                paramProduct = (feature.qualifiers["product"][0])
                
            if feature.type == "CDS":

                # length of product
                paramLength = round((feature.location.end - feature.location.start) / 3 - 1)

                tempRow = row.format(location=paramLocation,
                                     strand=paramStrand,
                                     length=paramLength,
                                     pid=pid+ int(i/2),
                                     gene=paramGene,
                                     synonym=paramSynonym,
                                     code="-",
                                     cog="-",
                                     product=paramProduct)
                pttRows.append(tempRow)

            elif feature.type in rntFeatureTypes:

                # length of product
                paramLength = round((feature.location.end - feature.location.start) - 1)

                tempRow = row.format(location=paramLocation,
                                     strand=paramStrand,
                                     length=paramLength,
                                     pid="-",
                                     gene=paramGene,
                                     synonym=paramSynonym,
                                     code="-",
                                     cog="-",
                                     product=paramProduct)
                rntRows.append(tempRow)

        # create files
        if(gbItem > 0):
            pttOutputFilePath = os.path.join(outputPath, "{}-{}.ptt".format(fileNameBase, gbItem))
            rntOutputFilePath = os.path.join(outputPath, "{}-{}.rnt".format(fileNameBase, gbItem))
        else:
            pttOutputFilePath = os.path.join(outputPath, "{}.ptt".format(fileNameBase))
            rntOutputFilePath = os.path.join(outputPath, "{}.rnt".format(fileNameBase))

        print("number of ptt rows:\t{}\t({})".format(len(pttRows), pttOutputFilePath))
        print("number of rnt rows:\t{}\t({})".format(len(rntRows), rntOutputFilePath))

        createFile(pttOutputFilePath, pttHeader, description, length, pttRows)
        createFile(rntOutputFilePath, rntHeader, description, length, rntRows)

        gbItem += 1

file = "GCF_005519465.1_ASM551946v1_genomic.gbff"
outputPath = ""
filePath = file
processFile(file, outputPath)
