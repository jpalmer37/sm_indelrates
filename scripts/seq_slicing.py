import os
from glob import glob
import re
import sys
from seqUtils import *


#GP120 Reference sequence file
with open(snakemake.input[1], 'r') as infile:
    gp120 = ''
    for line in infile:
        gp120 += line.strip("\n")

c_regions = [(0,390),  (588,885) , (993, 1152), (1254, 1377), (1410, 1532)]
#v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]
#modified
v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]

# builds an alignment index
def build_index(refseq):
    index = {}
    ri = 0 
    #Scan the reference for the variable region location
    for ai, x in enumerate(refseq):
        # ai = alignment index, literal position

        if x != '-':
            # otherwise alignment has a gap, do not increment reference index
            index.update({ri: ai})
            ri += 1

    return index

def get_regions(seq, index, seqtype='v'):

    positions = v_regions if seqtype == 'v' else c_regions

    slices = []
    for n1, n2 in positions:
        slices.append(seq[index[n1]:index[n2]])
    return slices


#Read and parse all subtype alignments
count = 0

filename = snakemake.input[0]

with open(filename) as handle:
    data = parse_fasta2(handle)

incorrect = []

#outputc = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/4_2_Conserved/" + filename + "_CR.fasta", "w")
with open(snakemake.output[0], "w") as cout, open(snakemake.output[1], 'w') as vout:

    for header, seq in data.items():
        #Extract the reference and query sequences
        ref, query = seq

        accno = header.split(snakemake.config['delimiter'])[snakemake.config['accno_pos']]

        # generate map from alignment to reference coordinates
         # reference index
        #{nucleotide #: actual position}
        index = build_index(ref)

        #FILTERING STEP
        #remove any sequences that contain a problematic variable region (more than 70% gaps)
        vout.write(accno + ",")
        cout.write(">" + header + "\n")

        vseqs = get_regions(query, index, 'v')
        
        gapcount = [float(seq.count("-")) / len(seq) for seq in vseqs]
        qcount = [float(seq.count("?")) / len(seq) for seq in vseqs]

        vseqs = [seq.replace("-","") for seq in vseqs]
        vseqs = ",".join(vseqs)

        vout.write(vseqs + "\n")

        cseqs = get_regions(query, index, 'c')

        cout.write(">"+header+"\n"+"".join(cseqs)+"\n")
        print(f"Finished slicing {header}")
# print(count)
# print(len(incorrect))
