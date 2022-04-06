import sys
from glob import glob
from seqUtils import *
import os


filename = snakemake.input[0] 
filename = os.path.basename(filename).split('_aligned')[0]


# First pass
# Find all positions with less than 95% gaps (whitelist)
with open(filename, "r") as fasta:

    transposed = transpose_fasta(convert_fasta(fasta))

    whitelist = []
    for pos, x in enumerate(transposed):

        gaps = x.count("-")

        freq = float(gaps)/len(x)

        if freq < 0.95:
            whitelist.append(pos)

# Second pass
# Write the positions in whitelist for all sequences
with open(filename, "r") as fasta, open(snakemake.output[0],'w') as outfile:

    data = parse_fasta3(fasta)

    for header, seq in data:
        fields = header.split(".")

        outseq = ''
        for n, char in enumerate(seq.copy()):

            if n in whitelist:
                outseq += char

        outfile.write('>' + header + "\n")
        outfile.write(outseq + "\n")