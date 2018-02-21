import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
import csv
import re
import scipy.io as sio
import scipy.stats
import random
import multiprocessing as mp
from math import log, exp, isnan
from copy import copy
from helpful_utils import reverse_complement
from collections import Counter
from interval import interval
from functools import partial
import argparse
import warnings
warnings.filterwarnings('ignore')


# Let's load in the Shendure model (Hexamer Additive Linear model, HAL), 
# determined an effect size for each exonic hexamer, as well as a model for the 
# exonic portion of the splice acceptor. The scoring functions and helper 
# functions are adapted from Rosenberg's github code used for the paper.

# DNA helper functions
bases = ['A','T','C','G']
def add_base(li):
    """Used in make_mer_list to add one more base to list"""
    new_li = []
    for s in li:
        for b in bases:
            new_li.append(s+b)
    return new_li
   

def make_mer_list(mer_len):
    """Makes a list of all n-mers"""
    li = bases
    for i in range(mer_len-1):
        li = add_base(li)
    return li


def logit(x):
    # refer to Wiki page on logit for asymptotic behavior at bounds
    if x == 0:
        return float('-Inf')
    elif x == 1:
        return float('Inf')
    else:
        return log(x) - log(1-x)


expit = lambda x: 1./(1.+exp(-x))


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


def score_exon_seq(seq, mer_scores, exonic_acceptor_scores, sd_scores, 
	mult_factor=1):
    if seq == '':
        return float('NaN')
    if 'N' in seq:
        return float('NaN')
    # need extra 5bp of right intron to properly score hexamers at exon/intron 
    # junction
    score = 0.
    # score first 3 bp of exon as part of splice acceptor
    score += exonic_acceptor_scores.ix[seq[:3]]*mult_factor
    # score rest of exon up until donor with exon hexamer scores
    for b in range(len(seq)-5-6-3): # don't score last 3bp
        score += mer_scores[seq[b:b+6]]*mult_factor
    # Score the hexamers overlapping with the exonic portion of the splice donor 
    # at -3, -2, -1:
    for b in range(3):
        score += sd_scores.ix[seq[len(seq)-8+b:len(seq)-8+6+b],b]*mult_factor
    return score


def make_exon_skipping_predictions(df, y_name, mult_factor=None):
    
    global exonic_mer6_scores
    global exonic_acceptor_scores
    global sd_scores

    # y_name is column name of splicing index to use for predictions
    if mult_factor is None:
        mult_factor = np.ones(len(df))*2.

    psi_pred_list = []
    for i in range(len(df)):
        # score mutant sequence
        mut_score = score_exon_seq(df.exon_seq.iloc[i], exonic_mer6_scores, 
          exonic_acceptor_scores, sd_scores, mult_factor[i])
        # score reference sequence
        ref_score = score_exon_seq(df.exon_seq_nat.iloc[i], exonic_mer6_scores, 
          exonic_acceptor_scores, sd_scores, mult_factor[i])
        psi_pred_list.append(expit(logit(df[y_name].iloc[i]) + mut_score - ref_score))
    
    df['PSI_pred'] = psi_pred_list
    df['DPSI_pred'] = df.PSI_pred - df[y_name]
    
    return df


def predict_with_error(scaling_factor, df):
    df = make_exon_skipping_predictions(df, y_name='v2_index', 
        mult_factor=np.ones(len(df))*scaling_factor)
    error = np.nansum((df.DPSI_pred-df.v2_dpsi)**2)
    return error


def optimize_scaling_factor(df):
    pool = mp.Pool(processes=4)
    params = np.arange(1,4,0.1)
    results = pool.map(partial(predict_with_error, df=df), params)
    scaling_factors = {params[i] : results[i] for i in range(len(params))}
    scaling_factors = pd.Series(scaling_factors)
    return scaling_factors.argmin()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Name of input file')
    parser.add_argument('outfile', help='Name of output file')

    args = parser.parse_args()
    # read in data
    exac_data = pd.read_table(args.infile, sep = '\t')
    # reformat so intron/exon length columns are int and not float, 
    # convert NA to 0 so we can use int
    exac_data[['start', 'end', 'intron1_len', 'exon_len', 'intron2_len']] = \
    	exac_data[['start', 'end', 'intron1_len', 'exon_len', 
    	'intron2_len']].fillna(0.0).astype(int)


    # load splice site models
    exonic_mer6_scores = pd.read_pickle('../../ref/Rosenberg_2015/exonic_mer6_scores.series')
    exonic_acceptor_scores = pd.read_pickle('../../ref/Rosenberg_2015/exonic_acceptor_scores.series')
    model = sio.loadmat('../../ref/Rosenberg_2015/model_full_data.mat')

    # grab model information for positions overlapping the splice donor (-3, -2, -1, +1)
    sd_scores = pd.DataFrame(index=make_mer_list(6),data=model['Mer_scores'][:4**6*8].reshape(4**6,8)[:,2:6])

    # extract exon sequence
    exac_data['exon_seq'] = exac_data.apply(lambda x : 
    	extract_exon_seq(x['strand'], x['intron1_len'], 
    		x['exon_len'], x['intron2_len'], x['original_seq'], 
    		extend=5, rc=False), axis=1)

    # extract exon sequence for natural sequence
    exac_data['exon_seq_nat'] = exac_data.apply(lambda x : 
        extract_exon_seq(x['strand'], x['intron1_len'], 
            x['exon_len'], x['intron2_len'], x['nat_seq'], 
            extend=5, rc=False), axis=1)

    # HAL only valid for exonic mutations, so let's subset our data and make 
    # our predictions.
    exac_exon_vars = exac_data[exac_data.label == 'exon']

    print "Determining optimal scaling factor through 10-fold cross-validation..."
    cross_val_inds = range(len(exac_exon_vars))
    random.seed(123)
    random.shuffle(cross_val_inds)
    cross_val_inds = np.array(cross_val_inds)
    # initialize
    cross_validated_scaling_factors = np.ones(len(exac_exon_vars))
    cross_val_group_size = len(exac_exon_vars)/10

    for i in range(10):
        print i
        test_inds = cross_val_inds[cross_val_group_size*i:cross_val_group_size*(i+1)]
        training_inds = np.array(list(set(range(len(exac_exon_vars)))-set(test_inds)))
        training_set = exac_exon_vars.iloc[training_inds]
        # assign factors learned on training set to the test set
        cross_validated_scaling_factors[test_inds] = optimize_scaling_factor(training_set)
        
    pickle.dump(cross_validated_scaling_factors, 
        open('../../processed_data/snv/HAL_cv_scaling_factors.pkl', 'wb'))

    # make predictions with the scaling factors for the held-out variants, so each
    # fold is assigned a scaling factor that was not trained on it
    exac_exon_vars = make_exon_skipping_predictions(exac_exon_vars, 'nat_v2_index',
        mult_factor=cross_validated_scaling_factors)

    exac_exon_vars.to_csv(args.outfile, sep='\t', na_rep='NA', index=False, 
    	quoting = csv.QUOTE_NONNUMERIC)
