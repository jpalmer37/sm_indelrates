from seqUtils import *
from glob import glob
import sys
import re 
import csv
import os 
import pandas as pd
from argparse import ArgumentParser

#listed sequences hinder the accurate formation of alignments
#delete the sequences in the given lists, and re-print the Conserved Region .fasta file

# parser = ArgumentParser()
# parser.add_argument('fasta',type=file)
# parser.add_argument('out',type=file)
# #parser.add_argument('out',type=lambda x: x if os.path.isdir(x) else NotADirectoryError)
# args = parser.parse_args()

# indir = snakemake.input[0]

# if not os.path.isdir(f'{args.out}/conserved/'):
#     os.mkdir(f'{args.out}/conserved/')
# if not os.path.isdir(f'{args.out}/variable/'):
#     os.mkdir(f'{args.out}/variable/')


blacklist =  {"01_AE":['KP411841'],
              "02_AG":['KP411843'],
              "C": ['KU319547','KP411838','MF373131','KU319550','KU319539','MF373138'],
              "B":['KP411824','KP411825','JQ403020','KT427845','DQ339453','KT427832']}

expr = re.compile("(.+)_CR.fasta")

subtype = os.path.basename(snakemake.input[0]).split("_CR")[0]
#csvpath = args.fasta.replace("conserved",'variable').replace("_CR.fasta",'_VR.csv')

try:
    fastafile = open(snakemake.input[0],"r")
    dataframe = pd.read_csv(snakemake.input[1],header=None)
except:
    print("ERROR: Reading fasta / csv files")
    sys.exit(1)
# Save the subtype name 


fasta = parse_fasta(fastafile)
# Iterate through the  
if subtype in blacklist.keys():
    for x in blacklist[subtype]:
        fasta.pop(x, None)

for accno in list(fasta.keys()):
    date = accno.split(snakemake.config['delimiter'])[snakemake.config['date_pos']]
    if date == snakemake.config['na_char']:
        print(f"Removing {accno}")
        fasta.pop(accno)

write_fasta(fasta, snakemake.output[0])


dataframe.columns = ['accno'] + ["V"+str(i) for i in range(1,6)]
dataframe.set_index("accno",inplace=True)


if subtype in blacklist:
    dataframe.drop(blacklist[subtype], axis=0,inplace=True,errors='ignore')

#dataframe.to_csv(f"{args.out}/variable/{subtype}_filter.csv")
dataframe.to_csv(snakemake.output[1])


