# extract splicing effects of SNVs from pre-computed file

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from helpful_utils import reverse_complement
import re
import scipy.io as sio
import scipy.stats
from math import log, exp, isnan
from copy import copy
import csv
from collections import Counter
from interval import interval
import warnings
warnings.filterwarnings('ignore')
import argparse


def extract_exon_seq(strand, intron1_len, exon_len, intron2_len, seq, extend=0, 
	rc=False):
    # extract exon sequence based on strand and intron/exon lengths. Grab n 
    # extra bp of downstream intron so hexamers at exon/intron border can be 
    # properly scored
    if not isinstance(seq, str):
        return ''
    
    seq = seq.upper()
    
    if strand not in ['+', '-', 1, -1, '1', '-1']:
        return ''
    if strand in ['-', -1, '-1']:
        if rc:
            seq = reverse_complement(seq)
        exon_seq = seq[intron2_len : intron2_len + exon_len + extend]
    else:
        exon_seq = seq[intron1_len : intron1_len + exon_len + extend]
    
    return exon_seq


def get_spanr(chr_num, position, mut_allele):
    # spidex is based on hg19
    coord = chr_num + ':' + str(position) + '-' + str(position)
    if not coord.startswith('chr'):
        coord = 'chr' + coord
    cmd = 'tabix ../../ref/spidex/spidex.tab.gz ' + coord + ' > tmp.tsv'
    os.system(cmd)
    output = pd.read_table('tmp.tsv', header=None, 
    	names=['chr', 'position', 'ref_allele', 'mut_allele', 'dpsi_max_tissue', 
    	'dpsi_zscore', 'gene', 'strand', 'transcript', 'exon_number', 'location', 
    	'cds_type', 'ss_dist', 'commonSNP_rs'])
    allele_match = output[output.mut_allele == mut_allele]
    
    return allele_match  


def predict_spanr(snps):
    no_spidex_match = []
    spanr_info = pd.DataFrame()
    for i in range(len(snps)):
        info = snps.iloc[i]
        snp_id = snps.id.iloc[i]
        spanr_output = get_spanr(info['chr'], info['snp_position'], 
        	info['alt_allele'])
        if len(spanr_output) == 0:
            no_spidex_match.append(snp_id)
        else:
            tmp = pd.DataFrame(spanr_output, 
            	columns=['chr', 'position', 'dpsi_max_tissue', 'dpsi_zscore', 
            	'location', 'ss_dist'])
            tmp['id'] = snp_id
            spanr_info = spanr_info.append(tmp)

    print "Number of variants with no SPIDEX match:", len(no_spidex_match)
    
    return spanr_info


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='Name of input file')
	parser.add_argument('outfile', help='Name of output file')

	args = parser.parse_args()

	exac_data = pd.read_table(args.infile, sep = '\t')
	# reformat so intron/exon length columns are int and not float, convert NA 
    # to 0 so we can use int
	exac_data[['start', 'end', 'intron1_len', 'exon_len', 'intron2_len']] = \
		exac_data[['start', 'end', 'intron1_len', 'exon_len', 
		'intron2_len']].fillna(0.0).astype(int)

	# subset to mutants only
	exac_snps = exac_data[exac_data.category == 'mutant']

	exac_spanr = predict_spanr(exac_snps)

	exac_spanr.to_csv(args.outfile, sep='\t', na_rep='NA', index=False)



