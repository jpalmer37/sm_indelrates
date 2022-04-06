from seqUtils import translate_nuc, parse_fasta
from gotoh2 import Aligner
from glob import glob
import os
import re
import multiprocessing as mp
import functools as ft

# class snake():
#     input = ['../data/input/A1.fasta','../data/reference/hxb2_gp120_sequence.txt']
#     output = ['../workflows/pair_aligned/A1_pairalign.fasta','../workflows/pair_aligned/A1_missing.txt' ]

# snakemake = snake()

#nucleotide version of the reference sequence
with open(snakemake.input[1], 'r') as infile:
    nt_ref = ''
    for line in infile:
        nt_ref += line.strip("\n")

#amino acid version of the reference sequence
aa_ref = translate_nuc(nt_ref,0)

#used for cutting off flanking regions
def get_boundaries(str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')

    # returns indices of without gap prefix and suffix
    res = [0, len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])

    if right:
        res[1] = len(str) - len(right[0])

    return res


stops = []
unequal = []

#nt alignment to remove any unwanted gene regions
nt_pair = Aligner()
nt_pair.set_model('HYPHY_NUC')
nt_pair.is_global = False
nt_pair.gap_open_penalty = 30
nt_pair.gap_extend_penalty = 10

#amino acid alignment
aa_pair = Aligner()
aa_pair.set_model('EmpHIV25')
aa_pair.gap_extend_penalty = 10
aa_pair.gap_open_penalty = 30
aa_pair.is_global = True

def pair_align(items):
    header, seq = items

    nt_result = nt_pair.align(nt_ref,seq)

    #finds the boundaries of the gene of interest
    left, right = get_boundaries(nt_result[0])
    
    #cuts the query sequence at the boundaries and changes it back to regular nt sequence
    nt_query = nt_result[1][left:right].replace('-','')

    #translate ntQry to amino acids
    aa_query = translate_nuc(nt_query,0)

    #skips any non functional sequences (you can check the genbank ID with the print line)
    if "*" in aa_query:
        return header, "NA", "NA"

    aa_result = aa_pair.align(aa_ref, aa_query)

    # convert the nucleotide sequences to lists so they can be edited (strings cant be mutated)
    final_ref = list(nt_ref)
    final_qry = list(nt_query)

    #reads through the amino acid alignment and adds codon gaps to the proper locations
    for i in range(len(aa_result[0])):
        if aa_result[0][i] == '-':
            final_ref[i*3:i*3] = ['-', '-', '-']

        if aa_result[1][i] == '-':
            final_qry[i*3:i*3] = ['-', '-', '-']


    if len(final_ref) != len(final_qry):
        #unequal.append(header)
        # print(nt_query)
        # print(nt_ref)
        return header, "NA", "NA"

    #aligned ref and query sequences
    final_ref = "".join(final_ref)
    final_qry = "".join(final_qry)

    return header, final_ref, final_qry

def multiprocess(data, processes, batch_size):
    with mp.Pool(processes) as pool:
        res = set()
        for val in pool.imap_unordered(pair_align, data.items()):
            print(val[0])
            
            res.update([val])

            if len(res) >= batch_size:
                print(f"Yielding {batch_size} alignments...")
                yield res
                res = set()

if __name__ == '__main__':
    with open(snakemake.input[0], 'r') as fasta_in:
        data = parse_fasta(fasta_in)

    with open(snakemake.output[0], 'w') as outfile, open(snakemake.output[1], 'w') as missing:
        for val in multiprocess(data, snakemake.threads, 50):
            for header, ref, qry in val:
                if ref == 'NA' or qry == "NA":
                    missing.write(header+"\n")
                else:
                    outfile.write(">" + header + '\n>ref\n' + ref + "\n>query\n" + qry + '\n')
    