'''
User modifiable global variables for the splicemod package

@author: dbgoodman
'''

import string
import os

from numpy import array

# import score

#===============================================================================
# general parameters
#===============================================================================

# where all this stuff should happen, usually one dir up from the src/ dir
topDir = os.getcwd()

# intron dropbox dir that contains sequence info
intronDataDir = os.path.join(topDir, 'data')

# not currently used
logFile = ('%s/parseFASTA.log' % topDir)

#####maxEnt Settings
# path to maxEnt executable
maxEntPath = os.path.join(topDir, "perl_utils/max_ent")

# donor input lengths for MaxEnt
maxEntBounds = {'me_donor': array([-3, 6]),
                'me_acceptor': array([-20, 3])}


# regions in acceptor and donor motifs that should not be mutated
maxEntInvariants = {'me_donor': [3, 4],
                    'me_acceptor': [18, 19]}

# what type of sequence is 5' of this motif type
maxEnt5Prime = {'me_donor': 'exon', 'me_acceptor': 'intron'}

# the filterScore determines at what score to print a cryptic motif
filterScore = 1

# these dicts will hold stored maxEnt scores for donors/acceptors
# this way, we won't have to recompute maxEnt scores with slow perl calls
maxEntScoreDicts = {'me_donor': {}, 'me_acceptor': {}}

######bs-finder Settings
findBSPath = topDir + "perl_utils"
findBSProg = "find_branchsite.pl"
# number of bases upstream of 3' acceptor to look for branch sites
bsUpstreamSearch = 50


######ppt-finder Settings
maxPPT3prime = 10  # maximum distance upstream of AG for ppt to end
minPPT3prime = 2  # minimum distance upstream of AG for ppt to end (incl AG)
minPPTpctCT = 0.5  # minimum percentage CT in PPT

######motif-finder Settings
motifDir = os.path.join(intronDataDir, 'motifs/')

######use pre-loaded pickled motifs
pickle_pfms = os.path.join(intronDataDir, 'motifs/motifs.pickle')

######context Settings
# NOTE: Currently not using...
context = ('ACCCGATTCAGCGAACCGCCTCGGTTTCCCTAACCCAATCCAGCCAGTAC',
           'TAATTTTCATAATTTGTTTTGTACTGAGTGCTGGCTAGTCAGATTACCTG')  # 2.2
context_size = 5


#===============================================================================
# mutate-motif Settings
#===============================================================================

max_mut_iter = 5  # number of times to mutate a single motif before stopping
mutate_record_iter = 5  # number of times to scan entire sequence iteratively
mutate_record_max = 5  # number of times to start over from original sequence

# when mutating motifs, if multiple motifs score below this value on maxent, then
# take one of them randomly. If all motifs are above this score, then just take
# the best one. Setting this to 0 means that if multiple mutations remove the
# original motif, then any of them are fine, but if none of them do, then only
# pick the best replacement
mutate_min_rand_choice = 1.5

correct_splice_signals = set()
# this will hold a unique set of splice signals that we want to keep


# see mutate.mutate_meta_feature() fxn for details on these
MUT_META_FINAL_PCT = .1
MUT_META_MAX_ITER = 3
MUT_META_MUT_PER_ITER = 2

# mammalian conservation settings
MAM_CONSERV_MIN_WINDOW_SCORE = 0.6

#####score_exons/generate_mutant Settings

# maxEnt score ranges to generate tuples in
gen_mut_ranges = [(-10, 0), (3, 5), (7, 8), (10, 20)]

# number of mutants to find per range (above)
gen_mut_count = 4

# this controls the amount of 'branching' when searching for mutants in a range,
# increasing will use more memory, but it might help find some mutant scores.
# 12 is good, lower will be faster, higher will get you more mutants if you are
# stuck. It corresponds to the number of single-nucleotide mutant branches to
# CREATE.
gen_mut_stored_sbps = 12

# similar implications to above, except it controls number of mutant branches to
# follow based on their closeness to the desired range after scoring.
gen_mut_follow_sbps = 5

######entropy Settings

# any string with a kolmogorov approx (computed with zlib) that is less than this
# number will not be kept when mutating motifs. This is to avoid long n-runs or
# short repeats. It's just an approximation of kolmogorov, of course, and doesn't
# work on very short strings.
use_entropy = True  # if False, don't calculate entropy changes

# currently not using improve_entropy
# if True, try to increase entropy if it starts too low
improve_entropy = True
# this it the minimum allowable entropy
kol_minscore = 0.5
# this is how much to expand around the motif when looking
kol_neighborhood = 4
# this is the size of the
kol_winsize = 8

#===============================================================================
# ensembl-exon Settings
#===============================================================================

# ensembl database settings

ens_release = 78

ens_lcl_db_dict = {'host':  "127.0.0.1",
                   'user':  "ensembl",
                   'passwd': "",
                   'port':   3306,
                   #'db':     "homo_sapiens_core_63_37"}
                   'db':     "homo_sapiens_core_78_38"}

ens_rmt_db_dict = {'host': 'ensembldb.ensembl.org',
                   'user': 'anonymous',
                   'port': 3306,
                   'db': "homo_sapiens_core_78_38",
                   'passwd': ''}

ens_data_dir = os.path.join(intronDataDir, 'ccds_ensembl/')
ens_gbk_dir = os.path.join(ens_data_dir, 'gbk/')
ens_fas_dir = os.path.join(ens_data_dir, 'fas/')

# this file is read from:
print ens_data_dir

ens_exon_fn = os.path.join(ens_data_dir, '78_38_CCDS_exons.all.txt')

# this file is written to:
ens_exon_stats = ens_data_dir + '1000_CCDS_exon_stats.txt'

cut_sites = ['GGCGCGCC',
             'TTAATTAA']

######common globals
programTemplate = string.Template("$path/$program")
# for convenience


#===============================================================================
# synthesis Settings
#===============================================================================
# Now we want to be able to synthesize more exons using longer oligos.

# the maximum exon/intron synthesis size, after removing primers
# (i.e. 170 if 200 base oligos, 190 if 230 base oligos)
chip_synth_length = 170

# total number of features to reserve for the controls, and to include in total
chip_reserved_features = 110
chip_feature_count = 55000 - chip_reserved_features

# minimum upstream and downstream intron length to include
chip_us_min = 40
chip_ds_min = 30

# this means the max exon size is chip_synth_length - chip_us_min - chip_ds_min:
# 200 - currently 90 bp or shorter
# 230 - currently 110 bp or shorter
# 300 - currently 180 bp or shorter

# there must be this much intronic space before/after exon
chip_intron_padding = 50

# number of exons to choose (randomly) to mutate
chip_number_to_mutate = 305







