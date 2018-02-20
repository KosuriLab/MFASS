'''
Created on Feb 17, 2011

@author: dbgoodman
'''

#external whole module imports
import sys
import cfg
import shutil
import re
import string
import collections
import Bio
import os
import warnings

import random

#local module imports
import util
import motif
import score
import mutate

#partial imports of classes
from copy import copy, deepcopy
from itertools import ifilter, chain, imap, combinations, product, repeat
from interval import interval
from operator import itemgetter
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from hashlib import md5


random_seed = 12345

def make_seq_feature(start, end, ftype, quals={}):
    '''
    create a sequence feature from a start, end, and a type. additionally you
    may include other fields, like note, label, evidence, citation, as a dict.
    '''

    seq_feature = SeqFeature(FeatureLocation(int(start), int(end)), strand= +1, type=ftype)
    seq_feature.qualifiers = quals
    seq_feature.qualifiers['source'] = ['splicemod']
    return seq_feature

#def extract_pos(obj):
#    '''
#    extract start/end tuples from a SequenceFeature object. Assumes that they
#    use ExactLocation as subobjects.
#    Can be called on a SeqFeature or a FeatureLocation (feature.location)
#    '''
#
#    #the function can be called either on a FeatureLocation or a SeqFeature, so
#    #test for type and handle both cases:
#
#
#    if isinstance(obj,Bio.SeqFeature.SeqFeature):
#        pos = obj.location
#    elif isinstance(obj,Bio.SeqFeature.FeatureLocation):
#        pos = obj
#    else:
#        raise TypeError('extract_pos() was given an object that was not ' + \
#                        'a SeqFeature or a FeatureLocation.')
#
#    return ( eval(str(pos.start)) , eval(str(pos.end)) )

def fix_sdb_genbank(seq_file):
    '''
    this little function ensures that there are 8 spaces padding the line
    numbers for genbank file. seqbuilder pads with 12, and stupid biopython
    crashes.
    '''

    shutil.move(seq_file, seq_file + "~")

    destination = open(seq_file, "w")
    source = open(seq_file + "~", "r")
    for line in source:

        line_match = re.match(r'^(\s+)(\d+)(.*)$', line)
        if line_match != None:
            destination.write(line_match.group(2).rjust(9, ' ') +
                              line_match.group(3) + '\n')
        else:
            destination.write(line)
    source.close()
    destination.close()

    return seq_file

#-------------------------------------------------------------------------------
#These classes will be added to Bio.SeqRecord via the util.magic_set decorator

@util.magic_set(SeqRecord)
def populate_attribs(self):
    '''
    This method adds everything to a new seqRecord that is needed for SPLICEMOD
    except the wiggle tracks, which must be added sequentially and so is done in
    the ensembl module.
    '''

    #this will remove all old splicemod features so we can find new ones
    self.remove_annotations()

    self.exon_list = []
    self.get_exon_list()
    self.correct_motifs = []
    self.find_all_motifs()
    self.get_correct_motifs()
    self.force_splice_signals()
    self.consolidate_all_motifs()
    self.snps_to_features()
    self.mut_sets = set()

@util.magic_set(SeqRecord)
def get_exon_list(self):
    '''
    find set of exonic intervals in a seq_record BUG: if exons aren't in
    order in gbk file, they wont be in exon_list, and this will crash the
    program.
    '''
    self.exon_list = []

    for sf in self.features:
        #if we encounter a CDS feature, change it to an exon
        if sf.type == 'cds':
            sf.type = 'exon'
            warnings.warn('Changed CDS at ' + str(sf.position) + ' to exon.\n')

        if sf.type == 'exon':
            self.exon_list.append(sf.location)

@util.magic_set(SeqRecord)
def add_conservation_features(self):
    '''
    will go through the wiggle track and make features of type conservation
    for intronic regions where conservation is higher than it should be for
    introns and is at least 6 bp wide
    '''

    tr = self.wigs['MamConserv']

    feature_locs = []
    idx = 0

    for scores in score.sliding_window(self.wigs['MamConserv'], 6, step=1):
        loc = (idx, idx + 6)
        if sum(scores) / 6 > cfg.MAM_CONSERV_MIN_WINDOW_SCORE:
            feature_locs.append(loc)
        idx += 1

    consol_flocs = interval.union([interval(*feature_locs)])

    #remove exons and splice signals:
    consol_flocs = consol_flocs & \
        interval(self.get_exon_and_signal_pos()).invert()

    for cfloc in consol_flocs:
        sta = int(cfloc[0])
        end = int(cfloc[1])
        if end - sta < 6: continue

        fscore = sum(tr[slice(sta, end)]) / (end - sta)
        self.features.append(make_seq_feature(sta, end, 'MamConserv',
                {'label'    :  ['MamConserv', ],
                 'evidence' : str(fscore)}))

@util.magic_set(SeqRecord)
def get_exon_and_signal_pos(self):
    feat_list = list(self.get_features('exon'))
    feat_list.extend(self.correct_motifs)

    ivl = interval(*[f.extract_pos() for f in feat_list])

    return (ivl[0][0], ivl[0][1])

@util.magic_set(SeqRecord)
def get_correct_motifs(self):
    '''
    go through every motif and set self.correct_motifs to the correct
    motif list
    '''

    self.correct_motifs = []

    for feat in self.features:
        if self.is_correct_motif(feat):
            self.correct_motifs.append(feat)

@util.magic_set(SeqRecord)
def replace_subseq(self, new_subseq, subseq_bounds):
    ''' This function edits a record given a short seq and its bounds in
        the sequence. It does not edit or add any features.

        Furthermore, the new_subseq string must be the same length as its
        bounds.
        In other words, if the original string was AAAAAA and the bounds
        were (2,5) then the new_subseq could be CCC and the updated record
        would contain the sequence AACCCA.
    '''
    if len(new_subseq) != (subseq_bounds[1] - subseq_bounds[0]):
        raise(ValueError('''This function cannot insert new sequence, only
                            replace existing sequence. The subseq_bounds
                            must be the same length as the new_subseq string
                         '''))

    seq_str = str(self.seq)
    new_seq = ''.join([seq_str[:subseq_bounds[0]],
                       new_subseq,
                       seq_str[subseq_bounds[1]:]])

    #change original sequence
    new_seq_obj = Bio.Seq.MutableSeq(new_seq, Bio.Alphabet.generic_dna)
    self.seq = new_seq_obj

@util.magic_set(SeqRecord)
def is_correct_motif(self, feat):
    '''
    this method looks through all exons in the exon list and check to see if
    feature should be there: if it is a donor, does it occur at the 3'
    border of an exon, and if it is an acceptor, does it occur at the 5'
    border.
    '''

    isCorrect = False

    # if it has already been marked,
    if ('function' in feat.qualifiers and
        feat.qualifiers['function'].count('correct splice signal')):
        return True

    # if it has not been marked, but is in the list
    if self.correct_motifs.count(feat):
        isCorrect = True


    #check if this is a splicemod splice signal
    else:

        #check to see if the 'label' qualifier value is one of the maxent
        #keys listed in the config file

        if  'label' not in feat.qualifiers:
            return False

        #get motif type
        if feat.type not in motif.motif_types:
            return False

        mt = motif.motif_types[feat.type]

        fiveprime = ''

        if 'fiveprime' in mt.attribs:
            fiveprime = mt.attribs['fiveprime']
        else:
            return False

        #based on the 5prime table, the sequence upstream should either be
        #an exon or an intron. if the latter, we will be checking the start
        #of the exon, i.e. the 5prime end, and the [0]th index of the tuple.
        if fiveprime == 'exon':
            exon_idx = 1 #end of exon, donor
        else:
            exon_idx = 0 #start of exon, acceptor

        motif_loc = feat.extract_pos()

        for exon in self.exon_list:
            ex_pos = exon.extract_pos()
            if motif_loc[0] == ex_pos[exon_idx] + mt.bounds[0]:
                isCorrect = True

            #in the case of an exon-exon boundary:
            if motif_loc[0] == ex_pos[not exon_idx] + mt.bounds[0]:
                isCorrect = False

    #finally, update the function field to reflect correctness
    if isCorrect:
        if 'function' not in feat.qualifiers:
            feat.qualifiers['function'] = []

        feat.qualifiers['function'].append('correct splice signal')

        return isCorrect

@util.magic_set(SeqRecord)
def remove_annotations(self):
    '''
    remove features generated by this program; they should all have source as
    'splicemod'
    '''

    new_features = []

    for feat in self.features:
        if not is_splicemod_motif(feat):
            new_features.append(feat)

    #sys.stderr.write("  Removed "+str(len(to_remove))+" old features.\n")
    self.features = new_features

@util.magic_set(SeqRecord)
def cut_out_splice_signals(self, feat):
    '''
    cut out correct splice signals from the interval and return a new interval
    (could have multiple components)
    '''

    spl_ivls = interval(*[f.extract_pos() for f in self.correct_motifs])

    if isinstance(feat, tuple):
        ivl = interval(feat)
    elif isinstance(feat, SeqFeature):
        ivl = interval(feat.extract_pos())

    return ivl & spl_ivls.invert()

@util.magic_set(SeqRecord)
def find_all_motifs(self):
    for motif_type in motif.motif_types.values():
        self.find_motifs(motif_type)

@util.magic_set(SeqRecord)
def find_motifs(self, seq_motif_type):
    '''associate successive n-mers with a motif_type object and add
    any found features to the record's feature list
    '''

    score_result_dict = \
        seq_motif_type.score(string.upper(str(self.seq)))


    nmers = score_result_dict['nmers']
    locations = score_result_dict['locations']
    scores = score_result_dict['scores']

    #names list is only given by some score types:
    if 'names' in score_result_dict:
        names = score_result_dict['names']
    else:
        names = []


    features = list()

    for i in range(len(locations)):

        #skip this site if it doesn't match the motif's filter criteria
        if not seq_motif_type.filter_score(scores[i]): continue

        #if there is a context attrib, skip if it's not in the right context
        if 'context' in seq_motif_type.attribs:

            exons = [interval(*[extract_pos(exon) for exon in self.exon_list])]
            exons = interval.union(exons)

            if seq_motif_type.attribs['context'] == 'exon':
                if interval(locations[i]) not in exons:
                    continue
            if seq_motif_type.attribs['context'] == 'donor_intron':
                donor_intron = interval([exons[-1][1], len(self.seq)])
                if interval(locations[i]) not in donor_intron:
                    continue

            if seq_motif_type.attribs['context'] == 'acceptor_intron':
                acceptor_intron = interval([0, exons[0][0]])
                if interval(locations[i]) not in acceptor_intron:
                    continue

        start = locations[i][0]
        end = locations[i][1]

        motif_bounds = seq_motif_type.bounds if \
            seq_motif_type.bounds else (0, len(nmers[i]))

        note_str = seq_motif_type.note_str(motif_bounds, nmers[i])

        feat = make_seq_feature(start, end, seq_motif_type.type,
                {'label'    :  [seq_motif_type.type, ],
                 'evidence' : str(scores[i]),
                 'note'     : [note_str, ],
                 'seq_motif_type': seq_motif_type.type}
                 )

        #add a motif name if this score type has them:
        if names:
            feat.qualifiers['instance_name'] = names[i]

        features.append(feat)

    self.features.extend(features)

@util.magic_set(SeqRecord)
def abridged_gb(self):
    '''if we are formatting 'gb', then don't print any features with a
       no_print motif attribute.'''

    abridged_feats = []

    for feat in self.features:
        if feat.type not in motif.motif_types:
            abridged_feats.append(feat)
        elif ('meta' in feat.qualifiers and feat.qualifiers['meta']) or \
          'no_print' not in motif.motif_types[feat.type].attribs:

            if 'me' in feat.type:
                feat.type = 'splice_site'

            abridged_feats.append(feat)

    new_record = deepcopy(self)
    new_record.features = abridged_feats

    gb_str = new_record.format('gb')
    del new_record
    return gb_str

@util.magic_set(SeqRecord)
def snps_to_features(self):

    for snp in self.annotations['snps']:

        if self.annotations['strand'] == 1:
            snp_start = snp.Location.Start - self.annotations['synth_start'] + 1
            snp_end = snp.Location.End - self.annotations['synth_start'] + 1
        else:
            snp_start = self.annotations['synth_end'] - snp.Location.End
            snp_end = self.annotations['synth_end'] - snp.Location.Start

        #skip if
        if snp_start < 0 or snp_end > len(self): continue

        feat = make_seq_feature(snp_start, snp_end, 'variation',
                {'label'    : 'variation',
                 'alleles'  : snp.Alleles,
                 'note'     : snp.Effect,
                 'pep'      : snp.PeptideAlleles,
                 'symbol'   : snp._get_symbol()})

        self.features.append(feat)

@util.magic_set(SeqRecord)
def consolidate_motifs(self, flist, mt):

    #return these lists:
    # stores new features that consolidate old ones that overlap
    consolidated_features = []
    # stores a list that tells which original features belong in which new one
    consolidated_list = {}


    ilist = [interval[extract_pos(feat.location)] for feat in flist]
    union = interval.union(ilist)

    #a dict of components and their constituent intervals
    cidict = {}

    #current interval
    civl = 0

    #for each component, check to see if this interval is in this connected
    #component. if it's not, go to the next component (since the ilist is
    #sorted when it is made).
    for comp in union.components:
        cidict[comp] = []
        while civl < len(ilist) and ilist[civl] in comp:
            cidict[comp].append(civl)
            civl += 1

    summed_scores = 0

    for comp in cidict:
        #sum up the components' scores
        fsum = 0
        fcount = 0
        flen = float(comp[0][1] - comp[0][0])
        for f_index in cidict[comp]:
            feat = flist[f_index]
            fsum += float(feat.qualifiers['evidence'])
            fcount += 1

        #calculate the feature mean, i.e. the average scores for all
        #sub-features, and the feature average, the average score
        #per base pair within the new metafeature
        fmean = float(fsum) / float(fcount) if fcount > 0 else 0
        flens = map(lambda i: i[0][1] - i[0][0], ilist)
        flen_avg = float(sum(flens)) / float(len(flens))
        favg = (fsum * flen_avg) / flen

        note_str = "len= %d; sum=%f; count=%d; mean=%f; avg=%f;" % \
            (flen, fsum, fcount, fmean, favg)

        #finally, make a new feature
        feat = make_seq_feature(comp[0][0], comp[0][1], mt.type,
            {'label'    :  [mt.type, ],
             'evidence' : str(favg),
             'note'     : [note_str, ],
             'seq_motif_type': mt.type,
             'nmers' : [str(self.seq[slice(*flist[fi].extract_pos())]) \
                           for fi in cidict[comp]],
             'meta'     : True})

        consolidated_features.append(feat)
        consolidated_list[feat] = []
        for f_index in cidict[comp]:
           consolidated_list[feat].append(flist[f_index])

    return (consolidated_features, consolidated_list)

@util.magic_set(SeqRecord)
def consolidate_all_motifs(self):
    '''for all motifs that can be consolidated, generate metafeatures that
       combine and average all of the contiguous scores.'''
    #if it's not made already, initialize a structure inside the seq_record obj
    # to keep track of the consolidated 'metafeatures' and their constituent
    # features
    if not hasattr(self, '_consolidated'):
        self._consolidated = {}

    #pick motif types that can be consolidated
    for mt in motif.motif_types.values():
        if 'can_consolidate' not in mt.attribs \
          or mt.attribs['can_consolidate'] is not True:
            continue

        self._consolidated[mt.type] = {}

        #list to hold the features
        flist = list(self.get_features(mt.type))

        (consolidated_features,
              consolidated_list) = self.consolidate_motifs(flist, mt)

        self.features.extend(consolidated_features)
        self._consolidated[mt.type] = consolidated_list


@util.magic_set(SeqRecord)
def generate_nmers(self, bounds):
    return util.generate_nmers(self.seq, bounds)

def add_missing_correct_splice_sites(self):
    '''
        force use of maxent on the flanks of the exon to get low splice
        scores
    '''

@util.magic_set(SeqRecord)
def add_wiggle_track(self, tr_name, wig_array):
    '''
    add a wiggle track to the wigs dict with track name as key to a SeqRecord
    '''
    if not hasattr(self, 'wigs'):
        self.wigs = {}
    self.wigs[tr_name] = wig_array

@util.magic_set(SeqRecord)
def __deepcopy__(self, memo_dict):
    new_record = copy(self)
    new_record.features = deepcopy(self.features, memo_dict)
    return new_record

@util.magic_set(SeqRecord)
def get_features(self, query, meta=False, aggregate=False, correct=None):
    '''
        returns an iterator on features given a certain query, a meta feature,
        a sequence motif type object, or a type string.
    '''

    #special cases for 'donor_intron' and 'acceptor_intron'
    if query == 'donor_intron':
        return ifilter(lambda f: len(self) in interval(f.extract_pos()),
                        self.get_features('intron'))
    if query == 'acceptor_intron':
        return ifilter(lambda f: 0 in interval(f.extract_pos()),
                        self.get_features('intron'))

    #if aggregate flag, then get all motif types that have that agg string
    if aggregate:
        query = ifilter(lambda m: 'aggregate' in m.attribs and \
                                   m.attribs['aggregate'] == query, \
                                 motif.motif_types.values())

    #if correct flag, then get all motif types that have the 'correct splice
    # signal' string
    if correct is True:
        is_correct_fxn = \
            lambda f: 'function' in f.qualifiers and \
                      'correct splice signal' in f.qualifiers['function']
        return ifilter(is_correct_fxn, \
                        self.get_features(query, \
                                          meta=meta, \
                                          aggregate=aggregate))
    if correct is False:
        not_correct_fxn = \
            lambda f: 'function' not in f.qualifiers or \
                      'correct splice signal' not in f.qualifiers['function']
        return ifilter(not_correct_fxn, \
                        self.get_features(query, \
                                          meta=meta, \
                                          aggregate=aggregate))

    sametype = lambda f, q: (f.type == q) or ('seq_motif_type' in f.qualifiers \
                            and f.qualifiers['seq_motif_type'] == q)

    if not isinstance(query, str) and isinstance(query, collections.Iterable):
        recurse_this = lambda q: self.get_features(q, meta, aggregate, correct)
        iter = chain.from_iterable(map(recurse_this, query))

    elif isinstance(query, str):
        iter = ifilter(lambda f: sametype(f, query), self.features)
    elif isinstance(query, motif.SeqMotifType):
        iter = ifilter(lambda f: sametype(f, query.type), self.features)
    elif isinstance(query, SeqFeature) and query.qualifiers['meta'] == True:
        iter = self._consolidated[f.qualifiers['seq_motif_type']][query]
    else:
        raise GetFeaturesQueryUnknownType

    return ifilter(lambda f: meta is f.is_meta(), iter)

@util.magic_set(SeqRecord)
def wig_for_feature(self, track, feat):
    '''get the wig scores for this feature's bases'''

    if feat not in self.features:
        raise FeatureNotFoundInSeqRecordException
    loc = feat.extract_pos()
    return self.wigs[track][slice(*loc)]

@util.magic_set(SeqRecord)
def wigs_for_ftype(self, ftype, track):
    return [self.wig_for_feature(track, feat)\
            for feat in self.get_features(ftype)]

@util.magic_set(SeqRecord)
def mutant_string(self, loc, tup):
    '''
    given a location and a set of mutant locations (loc/str tuples), give a
    mutant string of length location[1] - location[0] that make the changes
    specified.
    '''

    seq_str = str(self.seq)[slice(*loc)]
    return mutate.tups_to_str(seq_str, loc, tup)

@util.magic_set(SeqRecord)
def mutate_all_positions(self, loc):
    '''
    mutate every codon and/or nucleotide within feature bounds
    '''

    #first, make sure my self._mutant_locations dict is instantiated
    if not hasattr(self, '_mutant_locations'):
        self._mutant_locations = mutate.mutant_locations(self)

    #deal with interval() versus tuple inputs
    if isinstance(loc, interval):
        loc_ivl = loc
        loc = interval.hull([loc_ivl]).to_tuple()
    else:
        loc_ivl = interval(loc)

    ivl_len = loc_ivl.sum_len()

    (e_coords, i_coords) = mutate.get_motif_boundaries(loc, self)

    #mutant choices will be a list of random.choice lambda functions
    #that randomly chooses a different codon for every position or a different
    #nucleotide for every intronic base

    mutant_choices = set()

    for codon_loc in e_coords:

        #go through every codon in codon_loc
        for c_loc in range(codon_loc[0], codon_loc[1], 3):
            codon = str(self.seq[c_loc:(c_loc + 3)]).upper()
            #get other codons
            bckt = mutate.codon_back_table()
            fwdt = mutate.codon_fwd_table()
            other_codons = bckt[fwdt[codon]]
            other_codons = other_codons.difference((codon,))

            if len(other_codons) == 0:
                continue

            #convert these codons into mut tuples (cmut_tuples)
            # (one codon might be two or even three tuples)
            cmut_tuples = ()
            for other_cod in other_codons:
                cod_tup = ()
                for diff in util.str_diff(other_cod, codon):
                    diff_loc = (c_loc + diff, c_loc + diff + 1)
                    cod_tup += ((diff_loc, other_cod[diff]),)

                cmut_tuples += (cod_tup,)

            #finally store a lambda function that randomly chooses a
            #different codon for this position, using a unique-state
            #random generator
            rgen = random.Random()
            rgen.seed(random_seed ^ hash(cmut_tuples) ^ hash(loc))
            codon_choice = lambda cmt, rgen: lambda: rgen.choice(cmt)

            mutant_choices.add(codon_choice(cmut_tuples, rgen))

    for intron_loc in i_coords:

        intron_ivl = interval(intron_loc)
        mut_ivls = (interval(ml) for ml in self._mutant_locations.keys())
        loc_list = filter(itemgetter(1), \
                          [(ivl, ivl.overlaps(intron_ivl)) for ivl in mut_ivls])

        #change interval obj into loc tuple
        loc_tup = lambda loc: loc[0].to_tuple()
        #get the mutation set (the values) for a loc tuple
        loc_muts = lambda loc: self._mutant_locations[loc_tup(loc)]
        #expand the mutation set into individual mutations for a loc tuple

        loc_mset = lambda loc, rgen: \
            lambda: rgen.choice([((loc_tup(loc), i),) for i in loc_muts(loc)])

        #generate independently seeded random number gens for each pos
        rgens = [random.Random() for i in loc_list]
        [rg.seed((random_seed, loc)) for rg, loc in zip(rgens, loc_list)]

        pos_rnd_muts = map(lambda loc, rgen: loc_mset(loc, rgen), loc_list, rgens)

        mutant_choices.update(pos_rnd_muts)

    #now that we have a mutant choices list with one function for every
    # codon/nt, we need to create a generator that calls each function in the
    # list once only
    while True:
        yielded = set()
        next_mut = frozenset(
                        chain.from_iterable(map(lambda f: f(), mutant_choices)))
        seen_count = 0
        if next_mut not in yielded:
            yielded.add(next_mut)
            yield next_mut
        elif next_mut in yielded and seen_count < 20:
            seen_count += 1
        elif next_mut in yielded and seen_count >= 20:
            raise StopIteration

@util.magic_set(SeqRecord)
def force_splice_signals(self):
    '''
    if two splice signals are not recognized, then add them in despite their
    score
    '''
    mt_strs = ('me_splice_acceptor', 'me_splice_donor')
    mt_objs = [motif.motif_types[mt_str] for mt_str in mt_strs]

    for mt in mt_objs:
        if len(list(self.get_features(mt, correct=True))) > 0:
            continue
        fiveprime = mt.attribs['fiveprime']
        #based on the 5prime table, the sequence upstream should either be
        #an exon or an intron. if the latter, we will be checking the start
        #of the exon, i.e. the 5prime end, and the [0]th index of the tuple.
        if fiveprime == 'exon':
            exon_idx = 1 #end of exon, donor
        else:
            exon_idx = 0 #start of exon, acceptor

        ex_pos = self.get_features('exon').next().extract_pos()
        new_mtf_loc = tuple([i + ex_pos[exon_idx] for i in mt.bounds])
        new_mtf_str = self.seq[slice(*new_mtf_loc)].tostring()

        score = mt.score(string.upper(new_mtf_str))['scores'][0]

        note_str = mt.note_str(mt.bounds, new_mtf_str)

        feat = make_seq_feature(new_mtf_loc[0], new_mtf_loc[1], mt.type,
                {'label'    :  [mt.type, ],
                 'evidence' : str(score),
                 'note'     : [note_str, ],
                 'seq_motif_type': mt.type}

        )

        if 'function' not in feat.qualifiers:
            feat.qualifiers['function'] = []
        feat.qualifiers['function'].append('correct splice signal')
        self.features.append(feat)

@util.magic_set(SeqRecord)
def mutant_locations(self, loc, count=1, max=False):

    '''
    this function takes from the dict of all potential _mutation_locations a
    set of mutations that fall within a loc tuple. It returns an iterator that
    spits out potential mutants at these locations; it is randomized by position
    first, then by mutation.

    count is the number of mutations to return, setting to one returns all
    possible sequences off by one, setting to two returns all sequences with
    two mutations made, etc, etc.
    '''
    #first, make sure my self._mutant_locations dict is instantiatied
    if not hasattr(self, '_mutant_locations'):
        self._mutant_locations = mutate.mutant_locations(self)

    #deal with interval() versus tuple inputs
    if isinstance(loc, interval):
        if len(loc) == 0:
            return iter([])
        loc_ivl = loc
        loc = interval.hull([loc_ivl]).to_tuple()
    else:
        loc_ivl = interval(loc)

    #create an iterator that returns all keys for _mutant_locations that are
    #in this location range
    mut_ivls = (interval(ml) for ml in self._mutant_locations.keys())


    #now loc iter will output a non-random set of mutation locations which are
    #keys to the _mutation_locations dict
    loc_iter = ifilter(itemgetter(1), ((ivl, ivl.overlaps(loc_ivl)) for \
                                     ivl in mut_ivls))

    #change interval obj into loc tuple
    loc_tup = lambda loc: loc[0].to_tuple()
    #get the mutation set (the values) for a loc tuple
    loc_muts = lambda loc: self._mutant_locations[loc_tup(loc)]
    #expand the mutation set into individual mutations for a loc tuple
    loc_mset = lambda loc: ((loc_tup(loc), i) for i in loc_muts(loc))

    #put them all together for a randomized list of generators, one generator
    #for each loc tuple
    pos_mut_sets = map(lambda loc: (loc_mset(loc)), loc_iter)

    emit_sets = combinations(util.irandomize(
                    chain.from_iterable(pos_mut_sets), seed=random_seed), count)

    emit_sets = imap(frozenset, emit_sets)

    is_unique_pos = \
        lambda mset: (
            len([m[0] for m in mset]) == len(set([m[0] for m in mset]))
            and set(mset) not in self.mut_sets)

    mut_iter = util.irandomize(
                   ifilter(is_unique_pos,
                   util.irandomize(emit_sets,
                                   seed=random_seed)),
                   seed=random_seed)

    # if this feature overlaps exons
    #expand the motif to codons, so that we can check that mutants are
    # synonymous
    if interval(self.exon_list[0].extract_pos()).overlaps(loc_ivl):
        codon_loc = \
            (interval(mutate.expand_motif_to_codons(self, loc)) \
            & interval(self.exon_list[0].extract_pos())).to_tuple()

        #check all mutations for synonymousness
        seq_str = str(self.seq)[slice(*codon_loc)]
        is_synon = lambda seq_str, codon_loc: lambda mut_tups: \
           mutate.check_translation(\
               string.upper(mutate.tups_to_str(seq_str, codon_loc, mut_tups)),
               seq_str)

        is_synon = is_synon(seq_str, codon_loc)

        return util.irandomize(ifilter(lambda mut: is_synon(mut), mut_iter),
                               seed=random_seed)
    else:
         return mut_iter

@util.magic_set(SeqRecord)
def save_all_mutants(self, gb=True, fas=True):
    '''make a directory for this ens ID and save its original sequence and all
       mutants in two subdirs, FAS and GBK.
    '''

    paths_to_check = []
    seq_set = set()

    if gb:
        my_gbk_dir = cfg.ens_gbk_dir + self.id + '/'
        paths_to_check.extend([cfg.ens_gbk_dir, my_gbk_dir])

    if fas:
        my_fas_dir = cfg.ens_fas_dir + self.id + '/'
        paths_to_check.extend([cfg.ens_fas_dir, my_fas_dir])

    for path in paths_to_check:
        if not os.path.exists(path):
            os.mkdir(path)

    #add wild type to seq set
    seq_set.add(str.upper(self.seq.tostring()))

    #make one genbank file and one fasta file for each mutant

    if gb:
        seqrec_gbk = open(my_gbk_dir + self.id + '.gbk', 'w')
        seqrec_gbk.write(self.abridged_gb())
        seqrec_gbk.close()

    if fas:
        seqrec_fas = open(my_fas_dir + self.id + '.fas', 'w')
        seqrec_fas.write(self.format('fasta'))

    for mut_cat, mf_dict in self.mutants.items():
        for feat_name, mut_list in mf_dict.items():
            for mut in mut_list:

                mut_seq = str.upper(mut.seq().tostring())
                if mut_seq in seq_set:
                    warnings.warn('Duplicate sequence found! Skipping.')
                    continue
                else:
                    seq_set.add(mut_seq)

                if gb:
                    seqrec_gbk = open(my_gbk_dir + self.id
                                      + mut.fn() + '.gbk', 'w')
                    seqrec_gbk.write(mut.gb())
                    seqrec_gbk.close()

                if fas:
                    seqrec_fas.write(mut.fasta())
    if fas:
        seqrec_fas.close()

@util.magic_set(SeqRecord)
def generate_all_mutants(self):
    '''
        generate all mutants using mutant category
    '''
    #we want the items from the mut_cats to be populated first, then
    #we can make multi_muts from the pre-populated mut_cats
    all_mutant_generators = mutate.mut_cats.items()
    all_mutant_generators.extend(mutate.multi_mut_cats.items())

    #FOR EACH MUTANT CATEGORY...
    for mc_name, mut_cat in all_mutant_generators:
        print "{}".format(mc_name) + "=" * (80 - len(mc_name))

        #make a structure to store mutants of this mutant category
        self.mutants[mut_cat.cat] = {}

       #get a list of features to be mutated for this category
        mc_feats = list(mut_cat.get_mutant_fxn(self))

        #FOR EACH FEATURE IDENTIFIED BY/FOR THIS CATEGORY...
        for mcf in mc_feats:

            #if this mutant category has a skip_mutant function, and
            #this feature is true, then don't make this mutant. for
            #example, don't make a strong mutant if this feature is
            #already strong, don't make a weak mutant if this feature is
            #already weak, etc
            if  hasattr(mut_cat, 'skip_mut_fxn') \
              and mut_cat.skip_mut_fxn(mcf, self):

                name = mcf.name() if hasattr(mcf, 'name') else mc_name
                print "{} SKIPPED".format(name)\
                    + "-" * (72 - len(name))
                continue

            #set the mcf_name, None if this feature is not a SeqFeature object
            elif len(mc_feats) > 1 and mcf.name() is not None:
                print "{}".format(mcf.name())\
                    + "-"*(80 - len(mcf.name()))
                mcf_name = mcf.name()
            else:
                mcf_name = None

            #GENERATE MUTANTS FOR THIS FEATURE/CATEGORY COMBO
            fmuts = list(mut_cat.call_fxn(self, mcf))

            #update the set of mutant sets that we've already generated
            self.mut_sets.update([frozenset(fm.set) for fm in fmuts])

            #update the heirarchical structure of mutants in categories
            self.mutants[mut_cat.cat][mcf_name] = fmuts

            for num, mcm in enumerate(fmuts):
                mcm.num = num
                print "Mutant {}: {}".format(mcm.num, mcm)

                if mc_name is not 'variation' and not mcm.is_synonymous():
                    raise ValueError("MUTANT IS NOT SYNONYMOUS!")

#-------------------------------------------------------------------------------
#These classes will be added to Bio.SeqFeature via the util.magic_set decorator

@util.magic_set(SeqFeature)
def is_meta(self):
    return ('meta' in self.qualifiers and self.qualifiers['meta'])

@util.magic_set(SeqFeature)
def get_mtype(self):
    return motif.motif_types[self.type]

@util.magic_set(SeqFeature)
def is_splicemod_motif(self):
    '''
    Check that the current SeqFeature is a splicemod-annotated motif
    '''
    if not 'source' in self.qualifiers:
        return False

    elif (self.qualifiers['source'][0] == 'splicemod'):
        return True

    else: return False

@util.magic_set(SeqFeature)
def extract_score(self):
    '''returns the score associated with a feature'''

    if 'evidence' not in self.qualifiers:
        return None
    score = self.qualifiers['evidence']
    if type(score) == type([]):
        score = score[0]
    return float(score)

@util.magic_set(SeqFeature)
def __deepcopy__(self, memo_dict):
    new_feat = Bio.SeqFeature.SeqFeature(
                    Bio.SeqFeature.FeatureLocation(*self.extract_pos()),
                    type=str(self.type),
                    qualifiers=deepcopy(self.qualifiers, memo_dict),
                    strand=self.strand)
    return new_feat

@util.magic_set(SeqFeature)
def extract_pos(self):
    return (int(self.location.start.position), int(self.location.end.position))

@util.magic_set(SeqFeature)
def __len__(self):
    if not hasattr(self, '_len'):
        start, end = self.extract_pos()
        self._len = end - start
    return self._len

@util.magic_set(SeqFeature)
def __eq__(self, other):
    get_values = lambda feat: (feat.extract_pos(),
                               feat.type,
                               feat.extract_score())
    return get_values(self) == get_values(other)

@util.magic_set(SeqFeature)
def name(self):
    '''returns a pithy and hopefully unique name for this feature'''

    str = "{}_[{}:{}]".format(self.type, *self.extract_pos())
    if hasattr(self, 'num'):
        str += "-#{}".format(self.num)
    return str

@util.magic_set(FeatureLocation)
def extract_pos(self):
    return (int(self.start.position) , int(self.end.position))
