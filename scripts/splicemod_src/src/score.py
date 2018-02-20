'''
Created on Feb 17, 2011

@author: dbgoodman
'''

import subprocess
import glob
import string
import re
import entropy
#import motif

from numpy import array, where, bitwise_or, ones

import cfg
import feature

def call_find_branch_site(string):
    '''call find_branchsite and get a list of branch sites and consensus
    distances score of zero means consensus, every one NT away increases score
    output is "start   end   score   seq
    '''

    command = cfg.programTemplate.substitute(path=cfg.findBSPath,
                                             program=cfg.findBSProg)

    p = subprocess.Popen( ("perl",command,string),
                          stdin = subprocess.PIPE,
                          stdout = subprocess.PIPE,
                          close_fds = True,
                          cwd = cfg.findBSPath)
    stdout_text, stderr_text = p.communicate()

    #associate with fields
    bs_fields = ['start','end','score','seq']

    #split the output into line-by-line tuples and zip it into dicts
    zipsplit = lambda line: dict(
        zip(bs_fields,conv_ints_in_list(line.split())))
    return map(zipsplit,stdout_text.splitlines())

#def generate_nmers(self,*bounds):
#    '''generate a list of n-mers at each position in the seq_record, relative
#    to each position i, and remove any substrings that are too short
#    '''
#
#    #import pdb; pdb.set_trace()
#    #if one bound, start at 0 rel to i
#    if len(bounds) == 1:
#        bounds = array(0,*bounds)
#
#    #if two bounds, rel to i
#    elif len(bounds) == 2:
#        bounds = array(bounds)
#
#    #if nargs is not 1 or 2
#    else:
#        raise("relative bounds must be int or 2 ints")
#
#    seq_str = str(seq_record.seq)
#    str_lines = ''
#
#    for i in range(len(seq_str)):
#
#        #generate n-mers from bounds relative to position i
#        str_lines += "\n"+seq_str[slice(*(i + bounds))]
#
#    #remove short strings
#    str_list = clean_string_list(str_lines.splitlines(),sum(abs(bounds)))
#    return str_list


def find_ppt(seq):
    '''check for PPT sequence between a branch end site and a 5' acceptor
    must be a sequence w/ 50% CT
    '''

    seqa = array(list(seq))
    #find all valid start/end pairs (C/Ts)
    #valid start/ends
    start_end_ct = bitwise_or(seqa == 'T',seqa == 'C').nonzero()[0]

    best_coords = ()
    best_len = 0

    for start_i in range(len(start_end_ct)):

        start_l = start_end_ct[start_i]

        #maximum PPT 5' distance (after bs_end)
        #import pdb; pdb.set_trace()

        #if start_l <= bs_end:
        #    continue

        for end_i in range(start_i+1,len(start_end_ct)):

            end_l = start_end_ct[end_i]

            ppt_len = end_l - start_l + 1

            if ppt_len < best_len:
                continue

            dist_to_ag = (len(seq) - 1) - start_end_ct[end_i] - 2

            #minimum PPT 3' distance: -2, next to final AG
            #maximum PPT 3' distance -10, according to Kol 2005
            if not (cfg.maxPPT3prime - 2) >= dist_to_ag >= (cfg.minPPT3prime - 2):
                continue

            #perform pct_ct on all start/end pairs, take longest 'True' cell
            getPercentCT = lambda seq: (
                seq.count("C")+seq.count("T"))/float(len(seq))

            percentCT = getPercentCT(seq[start_l:(end_l+1)])

            if  percentCT > cfg.minPPTpctCT:
                best_len = ppt_len
                best_coords = (start_l,end_l)


    if best_len > 0:
        ppt = {'start': best_coords[0], 'end': best_coords[1]}
    else:
        ppt = {}
    return ppt

def find_strings(seq_record,str_sets):
    '''find string locations in multiple sets of motif strings
    string input format is
       str_sets[setName] = list of motif strings
    match output format is
       matched_sets[setName][stringName] = list of start/end tuples
    '''

    matched_sets = {}

    seq_str = str(seq_record.seq)
    for strSet in str_sets:
        matched_sets[strSet] = {}
        for string in str_sets[strSet]:
            matched_sets[strSet][string] = []

            lookStart = 0
            while True:
                strStart = seq_str.find(string,lookStart)
                #import pdb; pdb.set_trace()
                if strStart < 0:
                    break
                else:
                    locTuple = (strStart,strStart + len(string) - 1)
                    matched_sets[strSet][string].append(locTuple)
                    lookStart = strStart + 1

    return matched_sets

def load_motif_sets():
    '''
    create a motif set list from a list of motif txts in cfg.motifDir
    '''

    fnList = glob.glob(cfg.motifDir+"*.txt")
    motifSets = {}

    for fn in fnList:

        #name the set based on the file, stripping path and .txt
        setName = re.sub(r'.*/(.*)\.txt',r'\1',fn)

        #open the file and create a list of strings
        f = open(fn)
        motifSets[setName] = map(lambda x: x.upper(),f.read().splitlines())

    return motifSets

def annotate_splice_signals(seq_record):
    '''find all acceptors and donors with maxEnt'''

    #generate nmer lists for MaxEnt
    MEdonorList = generate_nmers(seq_record,*cfg.maxEntBounds['me_donor'])
    MEacceptorList = generate_nmers(seq_record,*cfg.maxEntBounds['me_acceptor'])

    donorFeatures = find_max_ent_motifs('me_donor',
                                     MEdonorList)
    seq_record.features.extend(donorFeatures)

    acceptorFeatures = find_max_ent_motifs('me_acceptor',
                                   MEacceptorList)

    seq_record.features.extend(acceptorFeatures)

    return seq_record

def conv_ints_in_list(vals):
    '''#convert string if it is an int, else leave it'''
    vals_out = []
    for val in vals:
        if (type(val) == type('') and val.isdigit()):
            vals_out.append(int(val))
        else:
            vals_out.append(val)
    return vals_out


def clean_string_list(strings,length):
    '''check that all strings in a list are of a certain length and clear out
    ones that are smaller
    '''
    strings = array(strings)
    strings[array([len(x) < length for x in strings])] = ''
    return strings

def score_seq_record(seq_record):
    '''
    sum all of the motif scores in a seq_record object as a rough way to gauge
    how well the mutation algorithms are doing.
    '''

    srScore = 0

    for feat in seq_record.features:
        if 'label' not in feat.qualifiers:
            continue
        if feat.qualifiers['label'][0].find('me_') == 0:

            if seq_record.correct_motifs.count(feat):
                continue

            fScore = max(0,feat.extract_score())
            srScore += fScore

    return srScore

def sliding_window(sequence,win_size,step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable.

    Courtesy of:
    scipher.wordpress.com/2010/12/02/simple-sliding-window-iterator-in-python/
    """

    # Verify the inputs
    try: it = iter(sequence)
    except TypeError:
        raise Exception("sequence must be iterable.")
    if not ((type(win_size) == type(0)) and (type(step) == type(0))):
        raise Exception("type(win_size) and type(step) must be int.")
    if step > win_size:
        raise Exception("step must not be larger than win_size.")
    if win_size > len(sequence):
        #raise Exception("win_size must not be larger than sequence length.")
        yield sequence

    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence)-win_size)/step)+1

    # Do the work
    for i in range(0,numOfChunks,step):
        yield sequence[i:i+win_size]

def entropy_score(seq):
    '''
    go over the sequence with a sliding window and return the window of size
    cfg.kol_winsize that has the lowest kolmogorov entropy score. if it is too
    low, this sequence shouldn't be used as a mutant (as it likely
    contains sequence of low complexity).
    '''

    min_entr = 1

    if seq < cfg.kol_winsize:
        raise(ValueError('The sequence given is too short to compute entropy!'))
    chunks = sliding_window(seq,cfg.kol_winsize)

    for chunk in chunks:
        chunk_entr = entropy.kolmogorov(chunk)
        if chunk_entr < min_entr:
            min_entr = chunk_entr

    return min_entr


