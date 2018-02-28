import pandas as pd
import numpy as np
import os
import pickle
import argparse
import sys
sys.path.append('./src')
import util
import feature
import motif
import Bio
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
import multiprocessing as mp
import csv
import math

pd.set_option('display.float_format', lambda x:'%f'%x)

# python scripts/score_lib.py exac/ref/exac_ref_formatted.txt splicemod/data/motifs/Ke2011/ ref/Rosenberg_2015/exonic_mer6_scores.series 4 exac/produced_data/exac_rescored.txt

def reverse_complement(seq):
    """
    Return the reverse complement of a nucleotide string
    """
    seq = seq.upper()
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    rc = ''.join([complement[nt] for nt in seq[::-1]])
    return rc


def extract_exon_seq(strand, intron1_len, exon_len, intron2_len, seq, extend=0, rc=False):
    # extract exon sequence based on strand and intron/exon lengths. Grab n extra bp of downstream intron so hexamers
    # at exon/intron border can be properly scored
    seq = seq.upper()
    
    if strand not in ['+', '-', 1, -1, '1', '-1']:
        return ''
    if rc:
        seq = reverse_complement(seq)
    if strand in ['-', -1, '-1']:
        exon_seq = seq[intron2_len : intron2_len + exon_len + extend]
    else:
        exon_seq = seq[intron1_len : intron1_len + exon_len + extend]
    
    return exon_seq


def score_exon(seq, score_dict, k):
    kmer_scores = [score_dict.get(seq[i:i+k], float('NaN')) for i in range(len(seq) - k + 1)]
    if len(kmer_scores) != 0:
        # remove nan from average
        score = np.mean([score for score in kmer_scores if not np.isnan(score)])
    else:
        return float('NaN')
    return score


def max_score_change(natural_seq, mutant_seq, score_dict, k, snv_position):
    if math.isnan(snv_position):
        return float('NaN')
    snv_position = int(snv_position - 1)
    # get scores for kmers that overlap SNV position
    mutant_kmers = [mutant_seq[i : i+k] for i in range(snv_position-k+1, snv_position + 1)]
    nat_kmers = [natural_seq[i : i+k] for i in range(snv_position-k+1, snv_position + 1)]
    mutant_scores = [score_dict.get(kmer, float('NaN')) for kmer in mutant_kmers]
    nat_scores = [score_dict.get(kmer, float('NaN')) for kmer in nat_kmers]
    # get difference in scores
    score_diff = [mutant_scores[i] - nat_scores[i] for i in range(len(mutant_scores))]
    # find biggest absolute change
    abs_scores = [abs(x) for x in score_diff]
    max_score_diff = score_diff[abs_scores.index(max(abs_scores))]
    return max_score_diff


def create_record(seq, intron1_len, intron2_len, exon_len, strand, name, rc):
    seq = seq.upper()
    if strand in [-1, '-1', '-']:
        if rc:
            seq = reverse_complement(seq)
        upstr_intron_len = intron2_len
        downstr_intron_len = intron1_len
    else:
        upstr_intron_len = intron1_len
        downstr_intron_len = intron2_len

    record = SeqRecord(Seq(seq), id=name)
    record.seq.alphabet = Bio.Alphabet.DNAAlphabet()
    record.annotations['upstr_intron_size'] = upstr_intron_len
    record.annotations['downstr_intron_size'] = downstr_intron_len
    record.annotations['exon_size'] = exon_len

    # create intron/exon features
    upstr_loc = (0, record.annotations["upstr_intron_size"])
    exon_loc = (upstr_loc[1], upstr_loc[1] + record.annotations["exon_size"])
    downstr_loc = (exon_loc[1], exon_loc[1] +
                   record.annotations["downstr_intron_size"])

    # upstream intron
    record.features.append(SeqFeature(FeatureLocation(*upstr_loc),
                                      type="intron", strand=strand))

    # exon
    record.features.append(SeqFeature(FeatureLocation(*exon_loc),
                                      type="exon", strand=strand))
    # downstream intron
    record.features.append(SeqFeature(FeatureLocation(*downstr_loc),
                                      type="intron", strand=strand))

    return record


def populate_attribs(record):
    # this exists as class function, define here so we can use it with multiprocessing
    record.exon_list = []
    record.get_exon_list()
    record.correct_motifs = []
    # sometimes find motifs doesn't work because the intron/exon is too short, 
    # like in the case of the imperfect variants. Catch these in a try clause
    try:
        record.find_all_motifs()
    except:
        return None
    # for some reason get_correct_motifs isn't working correctly, so just
    # manually do this
    # record.get_correct_motifs()
    for feat in record.features:
        if record.is_correct_motif(feat):
            record.correct_motifs.append(feat)

    try:
        record.force_splice_signals()
    except:
        return None
    try:
        record.consolidate_all_motifs()
    except:
        return None
    if 'snps' in record.annotations:
        record.snps_to_features()

    return record


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('lib_file', help='tab-separated library text file')
    parser.add_argument('--rev', help="optional, does library need to be reverse complemented", action='store_true')
    parser.add_argument('seq_name', help='column name of sequence', type=str)
    parser.add_argument('esr_folder', help='path to folder with ESE and ESS scores')
    parser.add_argument('mer6_file', help='exonic mer6 series')
    parser.add_argument('num_processes', help='Number of parallel processes', type=int)
    parser.add_argument('output_name', help='Name of output file')
    parser.add_argument('--debug', help='save records to .pkl for debugging', action='store_true')
    args = parser.parse_args()

    seq_name = args.seq_name
    data = pd.read_table(args.lib_file, sep='\t', dtype={'sub_id' : str})
    # reformat so intron/exon length columns are int and not float, convert NA to 0 so we can use int
    data[['intron1_len', 'exon_len', 'intron2_len']] = data[['intron1_len', 'exon_len', 'intron2_len']].fillna(0.0).astype(int)

    # # rename column
    # data.rename(index=str, columns={'seq' : 'sequence'}, inplace=True)

    data['exon_seq'] = data.apply(lambda x : extract_exon_seq(x['strand'], 
        x['intron1_len'], x['exon_len'], x['intron2_len'], x[seq_name], extend=0, rc=args.rev), axis=1)

    # Now that we have the exon sequences for our library, let's bring in the score 
    # matrices and give an effect size sum score for each sequence. A higher 
    # score means the hexamer is more likely to active nearby splice sites and 
    # promote exon inclusion (to be an ESE, exonic splicing enhancer). A negative 
    # score means the hexamer is an exonic splicing silencer and promotes exon skipping.
    print "Calculating average HAL exon score..."
    exonic_mer6_scores = pd.read_pickle(args.mer6_file)
    HAL_score_dict = exonic_mer6_scores.to_dict()
    data['avg_exon_effect_score'] = data.apply(lambda x : score_exon(x['exon_seq'], HAL_score_dict, 6), axis=1)
    
    # calculate max change in score between mutant and naturals
    # read in naturals
    # data['max_HAL_score_change'] = data.apply(lambda x: max_score_change(x['natural_seq'],
    #                                                                 x[seq_name], 
    #                                                                 HAL_score_dict, 6,
    #                                                                 x['rel_position']), axis=1)

    # calculate overall Ke score
    ESE_motifs = pd.read_table('./data/motifs/Ke2011/ESEseq.txt',
                          sep = '\t', header = None, names = ['seq', 'Ke_score'])
    ESS_motifs = pd.read_table('./data/motifs/Ke2011/ESSseq.txt',
                          sep = '\t', header = None, names = ['seq', 'Ke_score'])
    Ke_motifs = pd.concat([ESE_motifs, ESS_motifs])
    Ke_score_dict = Ke_motifs.set_index('seq').T.to_dict('records')[0]
    
    # # re-calculate exon sequence without extension
    # data['exon_seq'] = data.apply(lambda x : extract_exon_seq(x['strand'], 
    #     x['intron1_len'], x['exon_len'], x['intron2_len'], x['sequence'], extend=0, rc=args.rev), axis=1)
    
    data['Ke2011_avg_score'] = data.apply(lambda x : score_exon(x['exon_seq'], Ke_score_dict, 6), axis=1)

    # Let's format our library into `SeqRecord` objects from `biopython` to make 
    # the motif finding easier and compatible with existing code.

    data.strand = [-1 if strand == '-' else strand for strand in data.strand]
    data.strand = [1 if strand == '+' else strand for strand in data.strand]

    lib = data[(data.strand == -1) | (data.strand == +1)]
    print "Number of records:", len(lib)

    print "Creating records..."


    lib_records = lib.apply(lambda x : create_record(x[seq_name], x['intron1_len'], 
        x['intron2_len'], x['exon_len'], x['strand'], x['id'], rc=args.rev), axis=1)

    pool = mp.Pool(processes=args.num_processes)
    lib_records = pool.map(populate_attribs, lib_records)

    # filter out short records which returned None 
    lib_records = [record for record in lib_records if record is not None]
    print "Number of scored records: ", len(lib_records)

    # double check correct motifs are present
    for record in lib_records:
        if len(record.correct_motifs) < 2:
            for feat in record.features:
                if record.is_correct_motif(feat):
                    record.correct_motifs.append(feat)

    if args.debug:
        print "Dumping records..."
        pickle.dump(lib_records, open('lib_records.pkl', 'w'))
    
    # filter lib data frame to match scored records
    lib = lib[lib.id.isin([record.id for record in lib_records])]

    # let's calculate an average score for each feature type and store it as a column in the data frame (except splice
    # donor/acceptor, we'll deal with those separately)
    motifs = motif.motif_types
    for mtf in motifs:
        if not mtf.startswith('me'):
            feat_scores_all = [[feat.extract_score() for feat in list(record.get_features(motifs[mtf]))] for record in lib_records]
            lib[mtf+'_avg_score'] = [sum(feat_list)/len(feat_list) if len(feat_list) != 0 else None for feat_list in feat_scores_all]

    # merge lib with original data
    data = data.merge(lib[['id', 'Ke2011_ESS_avg_score', 'Ke2011_ESE_avg_score', 'Vlkr07_AICS_avg_score', 'Vlkr07_DICS_avg_score', 'RBPmats_avg_score']],
        how='outer', on='id')

    print len(data)

    # Let's create two columns, one each for the score of the correct splice acceptor and correct splice donor.
    data['correct_acc_seq'] = [''] * len(data)
    data['correct_acc_score'] = [float('NaN')] * len(data)
    data['correct_don_seq'] = [''] * len(data)
    data['correct_don_score'] = [float('NaN')] * len(data)
    data.set_index('id', drop=False, inplace=True)

    for record in lib_records:
        index = data[data.id == record.id].index.tolist()
        for feat in record.correct_motifs:
            if feat.type == 'me_splice_acceptor':
                data.set_value(index, 'correct_acc_seq', str(record.seq)[feat.extract_pos()[0]:feat.extract_pos()[1]])
                data.set_value(index, 'correct_acc_score', feat.extract_score())
            if feat.type == 'me_splice_donor':
                data.set_value(index, 'correct_don_seq', str(record.seq)[feat.extract_pos()[0]:feat.extract_pos()[1]])
                data.set_value(index, 'correct_don_score', feat.extract_score())

    data.to_csv(args.output_name, sep='\t', index=False, na_rep='NA', quoting=csv.QUOTE_NONNUMERIC)
