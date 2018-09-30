#!/usr/bin/python

'''
@author Kimberly Insigne
kiminsigne@gmail.com
Feb 2, 2017

This script takes as input a FASTQ file of reads and a reference FASTA file and
maps perfect DNA sequences to their matching reference sequence.

Output:
tab-separated text file, first column is Ensembl IDs and second column is 
corresponding list of sequences from the FASTQ file

'''

import re
import sys
import argparse
from itertools import islice
from helpful_utils import reverse_complement

def main():

    parser = argparse.ArgumentParser('home-made implementation of grep')
    parser.add_argument('-i', '--input', help='FASTQ reads file')
    parser.add_argument('-r', '--ref', help='FASTA reference file')
    parser.add_argument('-o', '--output', help='name of output file')
    args = parser.parse_args()

    reads = open(args.input, 'r')
    ref = open(args.ref, 'r')
    outfile = open(args.output, 'w')

    lib = {}

    while True:
        # parse reference file, reformat headers and convert to dictionary,
        # where key is sequence and value is header
        next_line = list(islice(ref, 2))
        if not next_line:
            break

        header = next_line[0].split()
        # mut = "nat_seq"             # default category is "natural sequence"
        # if "MD:Z" in next_line[0]:
        #     mut = "imperfect_seq"
        # if "cat=" in next_line[0]:     # parsing info only from splice-code ref files (chip 1)
        #     mut = header[6][:-1]    # parses mutant category information

        # take first 5 columns from header and "category info"
        # ID = header[0] + " " + header[1] + " " + header[2] + " " + header[3] + " " + header[4] + " " + mut   
        ID = header[0] + " " + header[1] + " " + header[2] + " " + header[3] + " " + header[4]  
        sequence = next_line[1].strip().upper() 
        lib[sequence] = ID
        # add reverse complement
        lib[reverse_complement(sequence)] = ID
        
    while True:
        # one entry is four lines in FASTQ format
        line = list(islice(reads, 4))
        if not line:
            break
        # sequence is second line of entry
        seq = line[1].strip()
        if seq in lib:
            outfile.write(lib[seq] + '\t' + seq + '\n')

    reads.close()
    ref.close()
    outfile.close()

main()
