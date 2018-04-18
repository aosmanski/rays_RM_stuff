#Import modules
import argparse
import pandas as pd
from pandas import DataFrame, read_csv
import sys
import os

#Set Arguments
def get_args():
	parser=argparse.ArgumentParser()
	parser.add_argument('-i', '--repeatmasker_output', type=str, 
						help='Repeat Masker output file', required=True)
	parser.add_argument('-l', '--location', type=str, help='location to start', 
						required=True)						
	parser.add_argument('-o', '--output', type=str, help='Output file name')
	args=parser.parse_args()
	RM=args.repeatmasker_output
	LOC=args.location
	OUT=args.output

	return RM, LOC, OUT

RM, LOC, OUT = get_args()
	
#Navigate to work folder.
os.chdir(LOC)

#read in the RepeatMasker output file. Include specified fields. 
df = pd.read_csv(RM, sep="\t", skiprows=(0,1,2), usecols=[4,5,6,8,9,10,15], names =["CONTIG", "START", "STOP", "ORIENT","TE_NAME","TE_FAM", "KIMURA"], header = None)

##############################################
### Format the dataframe to desired output ###

#Remove "kimura=" from the last field of the repeatmasker_output file.
df["KIMURA"] = df["KIMURA"].str.replace('kimura=', '')

#Find length of TE.
df['LENGTH'] = df['STOP'] - df['START'] - 1

#Split field 11 "TE_FAM" by the "/"... should create family and subfamily
df1 = pd.DataFrame(df.TE_FAM.str.split('/', 1).tolist(), columns = ["TE_FAMILY", "TE_SUBFAM"])

#Merge Dataframes
df2 = pd.concat([df, df1], axis=1)

#Change 'C' to '-' so bedtools can read dataframe
df2["ORIENT"] = df2["ORIENT"].str.replace('C', '-')

#Create new dataframe with the desired output fields
df3 = df2[["CONTIG", "START", "STOP", "ORIENT","TE_NAME","TE_FAMILY", "TE_SUBFAM", 'LENGTH', "KIMURA"]]

###############################################################
#### Creating element specific output files in .csv format ####

#Create outpute file of all SINEs
SINE=df3.loc[df3["TE_FAMILY"] == "SINE"]
SINE.to_csv(OUT + "_SINE.bed", sep="\t", header = False, index = False)

#Create output file of all DNA transposons
DNA=df3.loc[df3["TE_FAMILY"] == "DNA"]
DNA.to_csv(OUT + "_DNA.bed", sep="\t", header = False, index = False)

#Create outpute file of all rolling circle elements
DNA_RC=DNA.loc[DNA["TE_SUBFAM"] == "RC"]
DNA_RC.to_csv(OUT + "_DNA_RC.bed", sep="\t", header = False, index = False)

#Create outpute file of all non-LTR elements
NonLTR=df3.loc[df3["TE_FAMILY"] == "NonLTR"]
NonLTR.to_csv(OUT + "_NonLTR.bed", sep="\t", header = False, index = False)

#Create outpute file of all LTR elements
LTR=df3.loc[df3["TE_FAMILY"] == "LTR"]
LTR.to_csv(OUT + "_LTR.bed", sep="\t", header = False, index = False)

#Create outpute file of all Unknown elements
Unknown=df3.loc[df3["TE_FAMILY"] == "Unknown"]
Unknown.to_csv(OUT + "_Unknown.bed", sep="\t", header = False, index = False)
