from Bio import SeqIO
import math
import statistics
import argparse
import os
import pandas as pd
from pandas import *
import numpy


#To do list
#-1-Number of transcripts with BLAST hits in a reference proteome
#-2-Number of unique proteins in the reference proteome with a BLAST hit
#-3-Collapse factor (average number of transcripts that match the same protein sequence)


#Set Arguments
def get_args():
        parser=argparse.ArgumentParser()
        parser.add_argument('-i', '--blast', type=str,
                                                help='blast_output input file', required=True)
        parser.add_argument('-l', '--location', type=str, help='location to start',
                                                required=True)
        parser.add_argument('-o', '--output', type=str, help='Output file name')
        args=parser.parse_args()
        BLAST=args.blast
        LOC=args.location
        OUT=args.output

        return BLAST, LOC, OUT

BLAST, LOC, OUT = get_args()

#Navigate to work folder.
os.chdir(LOC)

#Create stats output file
STATS_OUTPUT=open(LOC + '/' + OUT + '_stats.txt', 'w')

############################################################################
####   Number of transcripts with BLAST hits in a reference proteome   #####
df = pd.read_csv(LOC + '/' + BLAST, sep="\t", usecols=[0,1], names =["Transcripts","Protein"], header = None)
NUM_HITS=len(df)
print("Total number of hits: " + str(NUM_HITS))

#Print to STATS_OUTPUT
STATS_OUTPUT.write("Total number of hits: " + str(NUM_HITS) + "\n")


###################################################
####   Number of unique proteins in the reference proteome with a BLAST hit   ####
UNIQ_PROT = set(df["Protein"])
print("Total number of unique proteins hit: " + str(len(UNIQ_PROT)))

#Print to STATS_OUTPUT
STATS_OUTPUT.write("Total number of unique proteins hit: " + str(len(UNIQ_PROT)) + "\n")


###################################################
### Collapsing factor ###
shite = df.groupby("Protein")["Transcripts"].count()
shite.tolist()
MEAN = shite.mean()
print("Average number of blast hits per protein: " + str(MEAN))

#Print to STATS_OUTPUT
STATS_OUTPUT.write("Average number of blast hits per protein: " + str(MEAN) + "\n")




#Do the collapsing factor using unigene


########################################################################
#######  Number of "unigenes" (unique genes in the transcriptome) #####
# The FASTA headers code for both unique genes and isoforms.
# For example: c9_g1_i1
#       c9_g1 = unique gene (unigene)
#       i1 = isoform
#
# Split the FASTA header and keep only the unique gene information

#Split field 11 "TE_FAM" by the "/"... should create family and subfamily
df1 = pd.DataFrame(df.Transcripts.str.split('_i', 0).tolist(), columns = ["unigene", "isoform"])

shite = df1.groupby("unigene")["isoform"].count()
shite.tolist()
MEAN = shite.mean()
print("Average number of blast hits per unigene: " + str(MEAN))

#Print to STATS_OUTPUT
STATS_OUTPUT.write("Average number of blast hits per unigene: " + str(MEAN) + "\n")
