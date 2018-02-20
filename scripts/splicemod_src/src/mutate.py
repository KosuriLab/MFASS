'''
Created on Feb 17, 2011

@author: dbgoodman
'''
from warnings import warn
from itertools import chain, islice, cycle, product, \
                      ifilter, combinations, repeat
from bisect import bisect, bisect_left
from numpy import array, delete
from copy import copy, deepcopy
from interval import interval
from operator import itemgetter
from blist import sortedset
from collections import OrderedDict

from Bio.Seq import Seq, MutableSeq
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation, SeqFeature

import cfg
import feature
import score
import util
import motif
import stats

import numpy
import types
import util
import sys
import random



numpy.seterr("raise")

#===============================================================================
# MUTANT CLASS
#===============================================================================
class Mutant:
    '''
    a class to hold information pertaining to a set of single nucleotide mutations
    that have a certain purpose including their category and the feature
    '''
    def __init__(self, sr, set,
                 score=None,
                 orig_score=None,
                 loc=None,
                 category=None,
                 notes={},
                 iter= -1):

        self.sr = sr
        self.set = set #the mutant set (really just a tuple of (loc,mut) tuples)
        self.orig_score = orig_score
        self.category = category #the mutant category object
        self.notes = notes
        self.iter = iter

        #loc can either be a feature or a location, in which case loc becomes
        #the extracted position and a new field, feature, is created for the
        #feature. Also create an ivl incase loc is disjoint
        if isinstance(loc, SeqFeature):
            self.loc = loc.extract_pos()
            self.ivl = interval(loc)
            self.feature = loc
        elif isinstance(loc, tuple):
            self.loc = loc
            self.ivl = interval(loc)
            self.feature = None
        elif isinstance(loc, interval):
            self.loc = interval.hull(loc)
            self.ivl = loc
            self.feature = None
        elif loc == None:
            self.loc = (0, len(sr))
            self.ivl = interval(0, len(sr))
            self.feature = None
        else:
            raise ValueError('Mutant __init__ failed, bad loc type!')

        #set the string
        self.str = sr.mutant_string(loc, set)

        #score this mutant if we can
        if score is None and category is not None:
            self.score = self.category.score_mutant_fxn(self)

        if orig_score is None \
                and self.feature is not None \
                and 'evidence' in self.feature.qualifiers:
            self.orig_score = float(self.feature.qualifiers['evidence'])

#from_set method is broken, does not cal __init__ correctly
#    @classmethod
#    def from_set(sr, set,
#                 loc= None,
#                 category= None,
#                 notes= {},
#                 iter= -1):
#        '''
#        this creates a mutant from a set. if no loc is given, the loc is the
#        entire seq record. (set = iterable of loc,nt tuples)
#        '''
#
#        if loc == None: loc = (0,len(sr))
#        return Mutant(sr, set, sr.mutant_string(loc,set), category, notes, iter)

    def with_new_mutants(self, new_set):

        #check that these new_set locations are not duplicates
        mut_dict = dict(self.set)
        if any([loc in mut_dict for loc, nt in new_set]):
            raise ValueError('Some mutant positions are already present!')

        new_set = frozenset(tuple(self.set) + tuple(new_set))

        return Mutant(self.sr,
                     new_set,
                     orig_score=self.orig_score,
                     loc=self.loc,
                     category=self.category)


    def __hash__(self):
        ''' mutant objects are hashable by their mutant tuple set'''
        return hash(self.set)

    def __key__(self):
        try:
            return (self.category.score_key_fxn(self.score, self.orig_score),
                    - len(self.set))
        except AttributeError:
            return self.__hash__()

    def __cmp__(self, other):
        return cmp(self.__key__(), other.__key__())

    def __lt__(self, other): return self.__key__() < other.__key__()
    def __le__(self, other): return self.__key__() <= other.__key__()
    def __eq__(self, other): return self.__hash__() == other.__hash__()
    def __ne__(self, other): return self.__hash__() != other.__hash__()
    def __gt__(self, other): return self.__key__() > other.__key__()
    def __ge__(self, other): return self.__key__() >= other.__key__()

    def __len__(self):
        return len(self.set)

    def __repr__(self):

        repr = "{}, {} changes\n".format(self.name(), len(self))
        repr += "\t{}\n".format(self.str)
        repr += "\t{}\n".format(self.loc_str())

        attribs_to_print = [(x, self.__dict__[x])
                            for x in ('score', 'orig_score')]
        fxns_to_print = [(x, getattr(self, x)())
                            for x in ('criteria_met',)]

        for k, v in chain(attribs_to_print, fxns_to_print, self.notes.items()):
            repr += "\t{}:\t{}\n".format(k, v)


        return repr

    def criteria_met(self):
        try:
            return self.category.keep_mutant_fxn(self.score, self.orig_score)
        except AttributeError:
            raise ValueError('''criteria_met() failed!
                                This mutant has no category,
                                so cannot meet criteria''')

    def name(self):
        '''a pithy and hopefully unique description of this mutant'''

        name_str = "{}_[{}:{}]".format(self.category.cat,
                                       self.loc[0],
                                       self.loc[1])

        #add a mutant number if we have one
        try:
            name_str += "-#{}".format(self.num)
        except AttributeError:
            pass
        finally:
            return name_str

    def fn(self):
        '''a pithy and hopefully unique and filename-asfe description of this
           mutant
        '''

        return "_{}_({}.{})-{}".format(self.category.cat,
                                       self.loc[0],
                                       self.loc[1],
                                       self.num)

    def stats(self):

        fields = [self.sr.id,
                  self.category.cat,
                  self.loc[0],
                  self.loc[1],
                  self.num,
                  self.str,
                  len(self),
                  '; '.join(['{}= {}'.format(*i) for i in self.notes.items()])]

        return '\t'.join(fields)

    @staticmethod
    def stats_header():

        fields = ['id', 'category', 'start', 'end', 'num', 'str', 'changes', 'notes']
        return '\t'.join(fields)

    def new_record(self):

        mut_record = copy(self.sr)

        #convert the seq into a mutable seq if it is not already
        if not isinstance(mut_record.seq, MutableSeq):
            mut_record.seq = mut_record.seq.tomutable()

        mut_record.seq[self.loc[0]:self.loc[1]] = self.str
        mut_record.populate_attribs()
        mut_record.features.append(self.tofeature())

        for feat in mut_record.features:
            if feat not in self.sr.features:
                feat.type = 'NEW_' + feat.type[0:11]

        if isinstance(self.feature, SeqFeature):
            old_feature = deepcopy(self.feature)
            old_feature.type = 'OLD_' + old_feature.type[0:11]
            mut_record.features.append(old_feature)

        mut_record.stats = stats.SRStats(mut_record)

        return mut_record

    def gb(self):

        return self.new_record().abridged_gb()

    def tofeature(self):

        qualifier_dict = {'label'       : self.category.cat,
                          'str'         : self.str,
                          'set'         : self.loc_str(),
                          'changes'     : len(self),
                          'orig_score'  : self.orig_score,
                          'new_score'   : self.score,
                          'criteria_met': self.criteria_met()}

        for k, v in self.notes.items():
            qualifier_dict[str(k)] = str(v)

        #ensure proper conversion to string
        for k, v in qualifier_dict.items():
            qualifier_dict[k] = str(v)

        feat = feature.make_seq_feature(self.loc[0], self.loc[1],
                                        'mutant', qualifier_dict)

        feat.num = self.num

        return feat

    def fasta(self):

        header_dict = OrderedDict([
                          ('cat'         , self.category.cat),
                          ('loc'         , ('{m.loc[0]}:{m.loc[1]}'
                                               ).format(m=self)),
                          ('str'         , self.str),
                          ('set'         , self.loc_str()),
                          ('changes'     , len(self)),
                          ('orig_score'  , self.orig_score),
                          ('new_score'   , self.score),
                          ('criteria_met', self.criteria_met())])

        for k, v in self.notes.items():
            header_dict[str(k)] = str(v)

        header_append = '; '.join(
            ['{}= {}'.format(str(k), str(v)) for k, v in header_dict.items()])

        mut_record = self.new_record()
        mut_record.description += ' ' + header_append
        fasta_str = mut_record.format('fasta')
        del mut_record
        return fasta_str

    def loc_str(self):
        return str([i for i in sorted(self.set)])

    def seq(self):
        '''Return the whole seqRecord with the mutation made'''

        mut_seq = self.sr.seq.tomutable()
        mut_seq[self.loc[0]:self.loc[1]] = self.str
        return mut_seq


    def is_synonymous(self):
        exon_loc = self.sr.get_features('exon').next().extract_pos()

        #if none of the mutations overlap the exon, then true
        if not interval(exon_loc).overlaps(
            interval.union([interval(s[0]) for s in self.set])):
               return True
        #otherwise, compare whole exon translations
        else:
            orig_str = self.sr.seq.tostring()[slice(*exon_loc)]
            mutant_str = \
               tups_to_str(orig_str, exon_loc, self.set)
            return check_translation(mutant_str, orig_str)

    def dist(self, other):
        '''
        computes the hamming distance between this mutant and another by
        taking the symmetric distance
        '''
        return set(self.set) ^ set(other.set)

    @staticmethod
    def separated_best(muts, best_size, min_separation):
        '''
        grab the best M mutants from a sortedset of mutants, but only take
        ones separated by at least N mutations
            pseudocode:
            1.  branch set starts with best mutant only
            2.  add addt'l muts starting from 2nd best on down only if none
                of the current muts in branch set are too close to it

            if min sep is one, just return top M from set, since they are unique
            by definition

            if there aren't enough separated mutants, add best unseparated ones

            if M = 0, keep taking the largest separated set

        '''
        if best_size == 0:
            return Mutant.separated(muts, min_separation)
        elif best_size >= len(muts):
            return muts

        sep_iter = Mutant.separated(muts, min_separation)

        best_set = sortedset(key=lambda mut: mut.__key__())
        best_set.add(sep_iter.next())

        while len(best_set) < best_size:
            try:
                best_set.add(sep_iter.next())
            #if we run out of the branch muts, then add back the best
            #unseparated mutants, but more separated mutants first
            except StopIteration:
                if min_separation > 2:
                    min_separation -= 1
                    sep_iter = ifilter(lambda mut: mut.criteria_met(),
                                       Mutant.separated(muts, min_separation))
                else:
                    sep_iter = Mutant.separated(muts, 1)

        return best_set

    @staticmethod
    def separated(muts, min_separation):
        '''
        returns an iterator based on muts (which is probably an iterator)
        that skips mutants that are not separated from previous mutants
        by at least min_separation changes
        '''

        def dist(mut1, mut2):
            set1, set2 = [mut.set if isinstance(mut, Mutant) else mut
                         for mut in (mut1, mut2)]
            set1, set2 = [set(set_i) if isinstance(set_i, tuple) else set_i
                         for set_i in (set1, set2)]
            #number of shared mutations should be at least min_separation
            #less than the smaller of the two
            return max(map(len, (set1, set2))) - len(set1 & set2)

        if min_separation <= 1:
            for mut in muts:
                yield mut

        yielded = set()

        def ok_to_add(mut, set):
            return (all([dist(mut, other) >= min_separation for other in set])
                    or len(set) == 0)

        for next_mut in muts:
            keep_this = ok_to_add(next_mut, yielded)
            if keep_this:
                yielded.add(next_mut)
                yield next_mut

#===============================================================================
# MUTANTCATEGORY CLASS AND CATEGORY INSTANCES
#===============================================================================
class MutantCategory:
    '''
    This will hold a 'type' or 'category' of mutation that will define things
    like how to score, which scores to keep, how many nucleotides to mutate,
    and how to pick mutants. It is basically a row in mutants_overview.xls,
    but all encoded into an object.

    pick_nt_fxn - the function used to return new mutant tuples
    min_msets_per_iter - minimum the number of new mutant sets to create per
                         iteration. If there arent enough unique mutants with
                         N changes, then the algorithm will lower the number
                         of mutants to generate the amount of diversity
                         necessary
    score_key_fxn - the function that acts as the sorting key for a mutant
    score_mutant_fxn - the function that scores a mutant
    keep_mutant_fxn - a t/f fxn that tells whether or not a mutant is 'good'
    add_muts_fxn - the function that will generate new mutant nts to add onto
                   an existing mutant
    new_per_branch - the maximum number of new mutants to add to an
                         existing mutant

    '''

    def __init__(self, **kwargs):

        mandatory_attrs = ['cat',
                           'get_mutant_fxn',
                           'motif_type',
                           'max_iter',
                           'nts_per_iter',
                           'muts_to_keep']

        default_attrs = {'pick_nt_fxn': MutantCategory.default_mut_loc,
                           'min_msets_per_iter': 1,
                           'score_key_fxn': lambda mc, new_sc, orig_sc:new_sc,
                           'score_mutant_fxn': lambda mc, mut: 0,
                           'keep_mutant_fxn': \
                                lambda mc, new_score, orig_score: True,
                           'add_muts_fxn': MutantCategory.default_add_muts,
                           'new_per_branch': 5,
                           'mut_iter_branch': 5,
                           'call_fxn': lambda self, sr, feat:
                                mutate_and_iterate(sr, feat, self),
                           'min_separation': 1
                          }

        missing_default = (lambda kw: lambda kv: kv[0] not in kw)(kwargs)
        all_args = chain(kwargs.items(),
                       ifilter(missing_default, default_attrs.items()))

        for k, v in all_args:
            if 'fxn' in k:
                setattr(self, k, types.MethodType(v, self, MutantCategory))
            else:
                setattr(self, k, v)

        for attr in mandatory_attrs:
            if not hasattr(self, attr):
                missing_attr = '"{}" attribute needed '.format(attr) \
                               + 'for a MutantClass object'
                raise TypeError(missing_attr)

    def default_add_muts(self, sr, mut, **kwargs):
        return add_muts_to_set(sr, mut, self, **kwargs)

    def meta_add_muts(self, sr, mut, **kwargs):
        '''
        create a complex disjoint interval to add mutant nts to, being the
        interval where clusters of motifs are found. if these motifs are splice
        signals, then remove real splice signals from mutation consideration

        shared by clst* and rmv*
        '''

        scoredict = motif.motif_types[self.motif_type].score(mut.str)

        if len(scoredict['scores']) == 0:
            #warn(("Mut set ({}) \n is devoid of motifs,"
            #      + " no new muts added.").format(mut.set))
            return []

        if 'splice' in self.cat:
            #remove real splice sites from our interval set
            for spl_loc in [f.extract_pos() for f in mut.sr.correct_motifs]:
                try:
                    spl_idx = scoredict['locations'].index(spl_loc)
                except:
                    pass
                else:
                    delete(scoredict['scores'], spl_idx)
                    delete(scoredict['locations'], spl_idx)


        scores = [float(s) for s in scoredict['scores']]
        locs = [(i + mut.loc[0], j + mut.loc[0]) for i, j in scoredict['locations']]

        if 'splice' in self.cat:
            locs, scores = zip(*[(loc, score) for loc, score in zip(locs, scores)
                         if score > 0])

            if len(locs) == 0:
                #this should be resolved somehow instead of raising
                raise ValueError("Only real splice signals are left")


        full_ivl = interval.union([interval(loc) for loc in locs])
        spl_ivls = interval(*[f.extract_pos() for f in sr.correct_motifs])
        mutable_ivl = full_ivl & spl_ivls.invert()

        if mutable_ivl.sum_len() == 0: return []
        return add_muts_to_set(sr, mut, self, loc=mutable_ivl, **kwargs)

    def meta_get_mutant(self, sr, mtype=None, meta=True):
        '''
        for the clust_* and rbpmat muts this function gets self.num_feat_to_mut
        features that will be mutated, picking the ones with the highest scores.
        '''
        if mtype == None:
            mtype = self.motif_type

        cl_number = self.num_feat_to_mut
        get_score = lambda f: float(f.qualifiers['evidence'])
        cmp_clusters = lambda fx, fy: cmp(get_score(fx), get_score(fy))

        cluster_list = list(sr.get_features(mtype, meta=meta))
        cluster_list.sort(cmp=cmp_clusters, reverse=True)
        return cluster_list[:min(len(cluster_list), cl_number)]

    def rmv_get_mutant(self, sr, mtype=None, meta=True):
        '''
        for the clust_* and rbpmat muts this function gets self.num_feat_to_mut
        features that will be mutated, picking the ones with the highest scores.
        '''

        #dont return any mutants if there are no extra splice sites
        if 'splice' in self.cat:
            orig_score = sr.stats['motif'][str(self.motif_type)]['sum']
            good_score = sum(
                [f.extract_score() for f in sr.get_features(
                                                self.motif_type,
                                                correct=True)])
            if orig_score == good_score:
                return []

        if 'context' in motif.motif_types[self.motif_type].attribs:
            context = motif.motif_types[self.motif_type].attribs['context']
            context_feats = list(sr.get_features(context))
            if len(context_feats) > 1:
                raise ValueError('Too many context features!')
            if len(context_feats) == 0:
                raise ValueError('No context features!')
            return [context_feats[0].extract_pos()]
        else:
            return [(0, len(sr))]


    def rmv_mut_loc(self, sr, feat, count=1):
        '''
        for the initial mutant set, add mutants to an empty mutant!
        '''

        orig_score = (sr.stats['motif'][str(self.motif_type)]['sum']
                      - sum([f.extract_score() for f in sr.get_features(
                                                            self.motif_type,
                                                            correct=True)]))

        if 'splice' in self.cat:
            good_scores = sum(
                [f.extract_score() for f in sr.correct_motifs])
            orig_score = orig_score - good_scores


        empty_mut = Mutant(sr, frozenset(),
                           loc=feat,
                           orig_score=orig_score,
                           category=self)

        muts = self.add_muts_fxn(sr, empty_mut,
                                 nts_to_add=self.nts_per_iter,
                                 num_new_sets=self.new_per_branch)

        return (mut.set for mut in muts)

    def rmv_call_fxn(self, sr, feat):
        orig_score = (sr.stats['motif'][str(self.motif_type)]['sum']
                      - sum([f.extract_score() for f in sr.get_features(
                                                            self.motif_type,
                                                            correct=True)]))
        return mutate_and_iterate(sr, feat, self, orig_score=orig_score)

    def rmv_score_mutant(self, mut):
        '''
        same as meta_score_fxn but we are taking sums not fsums
        '''
        scoredict = motif.motif_types[self.motif_type].score(mut.str)
        if len(scoredict['scores']) == 0:
            return 0

        scores = [float(s) for s in scoredict['scores']]

        if 'splice' in self.cat:
            #remove real splice sites
            for spl_loc in [f.extract_pos() for f in mut.sr.correct_motifs]:
                try:
                    spl_idx = scoredict['locations'].index(spl_loc)
                except:
                    pass
                else:
                    scores.pop(spl_idx)

            #only take good scores
            fsum = sum(filter(lambda s: s > 0, scores))
            return fsum
        else:
            fsum = sum(scores)
            return fsum

    def meta_score_mutant(self, mut):

        scoredict = motif.motif_types[self.motif_type].score(mut.str)

        if len(scoredict['scores']) == 0:
            return 0

        scores = [float(s) for s in scoredict['scores']]

        if 'splice' in self.cat:
            #remove real splice sites
            for spl_loc in [f.extract_pos() for f in mut.sr.correct_motifs]:
                try:
                    spl_idx = scoredict['locations'].index(spl_loc)
                except:
                    pass
                else:
                    scores.pop(spl_idx)

            #only take good scores
            fsum = sum(filter(lambda s: s > 1, scores))
            return fsum
        else:
            fsum = sum(scores)
            flen = len(mut.str)
            #sum of lengths of individual nmers
            flens = map(lambda i: i[1] - i[0], scoredict['locations'])

            flen_avg = float(sum(flens)) / float(len(flens))
            favg = (fsum * flen_avg) / flen
            return favg

    def conserved_get_mutant(self, sr):
        '''
        gets conserved features, up to self.num_feat_to_mut, based on total
        conservation (sum), which is a balance between choosing conserved
        features based on size of feature and amount of conservation
        '''

        consv_feats = list(sr.get_features('MamConserv'))

        #this is the key function, sorting features by the sum of their
        sort_consv_sum = (lambda sr: lambda feat: \
            sum(sr.wigs['MamConserv'][slice(*feat.extract_pos())]))(sr)

        consv_feats.sort(key=sort_consv_sum, reverse=True)
        return consv_feats[:min(len(consv_feats), self.num_feat_to_mut)]

    def conserved_mut_loc(self, sr, feat, count=1):
        '''
        for the conserved mutant clusters, this picks the best
        $count bases & calls sr.mutant_locations to mutate them.
        '''
        #get the values for the wig track at this conservation feature loc
        loc = feat.extract_pos()
        values = sr.wigs['MamConserv'][slice(*loc)]
        #associate them with their locations
        val_tups = zip(range(*loc), values)
        #sort them by value
        best_tups = sorted(val_tups, key=itemgetter(1), reverse=True)

        mut_tuples = ()

        #loop through best tups self.nts_per_iter times and return a set of
        #random 1bp mutants
        for i in range(count):
            tup_loc = (best_tups[i][0], best_tups[i][0] + 1)
            mut_tuples += tuple(sr.mutant_locations(tup_loc).next())

        return (frozenset(mut_tuples),)

    def default_mut_loc(self, sr, feat, count=1, mut_iter=None):
        '''
        this is the default function for mutant categories to get mutants, it
        is just a wrapper for sr.mutant_locations
        '''

        return Mutant.separated(
                   sr.mutant_locations(feat.extract_pos(), count),
                   self.min_separation)

    def random_mut_loc(self, sr, feat, count=1):
        '''
        this is the default function for mutant categories to get mutants, it
        is just a wrapper for sr.mutant_locations
        '''
        rands = Mutant.separated(
                    sr.mutant_locations(feat.extract_pos(), count),
                    self.min_separation)

        return [rands.next() for i in range(self.muts_to_keep)]

    def splice_mut_loc(self, sr, feat, count=1, invert=False, full=False):
        '''
        this function takes a splice signal feature and gets mutant nucleotides
        in a disjoint interval that excludes conserved motif locations, as
        specified by the motif type object associated with the feature. if
        invert is set, then it ONLY returns mutants that are conserved motif
        locations.
        '''

        #get mtype
        mtype = motif.motif_types[feat.qualifiers['seq_motif_type']]
        #get interval
        ivl = interval(feat.extract_pos())
        #get invariant intervals
        start = ivl.extrema[0][0]
        make_ivl_tup = lambda iv: (iv + start, iv + start + 1)
        invar = interval(*map(make_ivl_tup, mtype.invariant))

        #invert and difference the invariant from the splice signal interval
        if not full:
            if not invert:
                final_ivl = ivl & invar.invert()
                #if invert is on, ONLY take conserved nucleotides
            else:
                final_ivl = invar
                count = 2
        else:
            final_ivl = ivl

        return Mutant.separated(sr.mutant_locations(final_ivl, count),
                                self.min_separation)

    def rbp_score_mutant(self, mut):
            scoredict = motif.motif_types['RBPmats'].score(mut.str);
            if len(scoredict['locations']) == 0 \
              or scoredict['locations'][0] == 0:
                return 0
            else:
                return scoredict['scores'][0]

    def snp_get_mutant(self, sr):
        #filter out snps that have bad effects
        snps = ifilter(lambda snp: not stats.has_effect(stats.bad_effects, snp),
                       sr.get_features('variation'))

        #filter out snps that change the length of the sequence

        def no_length_change(snp):
            nts = ('A', 'G', 'C', 'T')
            pos_len = snp.extract_pos()[1] - snp.extract_pos()[0]
            alleles = snp.qualifiers['alleles'].split('/')
            #are the lengths of all alleles the same?
            if not len(set([len(a) for a in alleles])) == 1:
                return False
            #do all the alleles consist of ATGC? (i.e. no '-')
            elif not all([allele_nt in nts \
                          for allele_nt in chain.from_iterable(alleles)]):
                return False
            elif not all([len(allele) == pos_len for allele in alleles]):
                return False
            else:
                return True

        return ifilter(no_length_change, snps)

    def snp_mut_loc(self, sr, feat, count=1):
        strand = sr.annotations['strand']

        mut_tuples = set()

        #find wild type allele at this locus, doing a rev-trans if necessary
        wt_allele = sr.seq[slice(*feat.extract_pos())]
        if strand == -1:
            wt_allele = wt_allele.reverse_complement()
        wt_allele = wt_allele.tostring()

        #get all alleles from the snp allele string
        alleles = set(feat.qualifiers['alleles'].split('/'))

        #remove the wild type
        if not wt_allele in alleles:
            raise ValueError(
                'Wild type allele {} not in allele set {}'.format(wt_allele,
                                                                  alleles))
        mut_alleles = alleles - set(wt_allele)

        if strand == -1:
            mut_alleles = set([Seq(ma).reverse_complement().tostring()
                               for ma in mut_alleles])

        start, end = feat.extract_pos()

        for mut in mut_alleles:
            mut_set = ()
            for i, mut_nt in enumerate(mut):
                pos = (start + i, start + i + 1)
                mut_set += ((pos, mut),)
            mut_tuples.add(mut_set)

        return frozenset(mut_tuples)

    def csplice_get_mutant(self, sr):

        not_correct = lambda f: \
            'function' not in f.qualifiers or \
            'correct splice signal' not in f.qualifiers['function']
        all_csplice = list(ifilter(not_correct, sr.get_features(\
                        self.motif_type)))

        get_score = lambda f: float(f.qualifiers['evidence'])
        cmp_splice = lambda fx, fy: cmp(get_score(fx), get_score(fy))

        all_csplice.sort(cmp=cmp_splice)
        return all_csplice[:min(len(all_csplice), self.num_feat_to_mut)]

    def aggr_get_mutants(self, sr):

        if self.motif_type is 'both':
            return (interval((0, len(sr))).to_tuple(),)
        else:
            return sr.get_features(self.motif_type)

    def aggr_mut_loc(self, sr, feat, count):
        '''call mutate_all_positions for the region
           a mutation to each mutable site, (except splice signals)
        '''

        #cut out splice signals from the interval
        final_ivl = sr.cut_out_splice_signals(feat)
        #send potentially disjoint interval to mutate_all_positions()
        return Mutant.separated(sr.mutate_all_positions(final_ivl),
                                self.min_separation)

    def paired_mut_loc(self, sr, feat, count):
        ''' this merges multiple mutants into one, with a new mutant category'''

        if isinstance(self.motif_type, tuple):
            mut_pairs = zip(sr.mutants[self.motif_type[0]][None],
                            sr.mutants[self.motif_type[1]][None])
        else:
            mut_pairs = zip(*(sr.mutants[self.motif_type].values()))

        for mp in mut_pairs:
            #check that no mutant sets are overlapping
            comb_iter = combinations([m.loc for m in mp], 2)
            if any([interval(i) in interval(j) for i, j in comb_iter]):
                raise ValueError('These intervals overlap!')

            yield frozenset(chain.from_iterable([m.set for m in mp]))


#------------------------------------------------------------------------------
mut_cats = OrderedDict() #a dict of mutant categories


multi_mut_cats = {} # a dict of mutant categories that consist of combinations
                    # of mutants in the mut_cats dict above.

#------------------------------------------------------------------------------
#MUTANT CATEGORIES - they are put in the mut_cats dict and are called with
#                    mutate_and_iterate()

#===============================================================================
# SNPs - SNP mutants
#===============================================================================
mut_cats['variation'] = \
    MutantCategory(\
        cat='variation',
        score_mutant_fxn=lambda mc, mut: 0,
        get_mutant_fxn=MutantCategory.snp_get_mutant,
        pick_nt_fxn=MutantCategory.snp_mut_loc,
        motif_type='variation',
        max_iter=1,
        nts_per_iter=1,
        muts_to_keep=1,
        mut_iter_branch=0,
        num_feat_to_mut=10)

#these categories aggregate mutants made from the mut_cats, dict, so we have to
#generate them after, by appending their keys in the ensembl module
#===============================================================================
# clust_* - METAMUTANT CLUSTERS
#===============================================================================
for meta_motifs in ['Ke2011_ESE', 'Ke2011_ESS', 'Vlkr07_AICS', 'Vlkr07_DICS']:
    mut_cats['clst_' + meta_motifs] = \
        MutantCategory(\
            cat='clst_' + meta_motifs,
            keep_mutant_fxn=lambda mc, new_score, orig_score: \
                abs(orig_score) * 0.1 > abs(new_score),
            score_mutant_fxn=MutantCategory.meta_score_mutant,
            get_mutant_fxn=MutantCategory.meta_get_mutant,
            add_muts_fxn=MutantCategory.meta_add_muts,
            motif_type=meta_motifs,
            max_iter=20,
            nts_per_iter=1,
            muts_to_keep=2,
            mut_iter_branch=10,
            num_feat_to_mut=2,
            min_msets_per_iter=2,
            score_key_fxn=lambda mc, new_score, orig_score: abs(new_score),
            min_separation=2)

#===============================================================================
# rmv_* - REMOVE ALL METAMUTANT CLUSTERS
#===============================================================================
for meta_motifs in ['Ke2011_ESE', 'Ke2011_ESS', 'Vlkr07_AICS', 'Vlkr07_DICS',
                    'me_splice_donor', 'me_splice_acceptor']:
    mut_cats['rmv_' + meta_motifs] = \
        MutantCategory(\
            cat='rmv_' + meta_motifs,
            pick_nt_fxn=MutantCategory.rmv_mut_loc,
            keep_mutant_fxn=lambda mc, new_score, orig_score: \
                abs(orig_score) * 0.1 > abs(new_score),
            score_mutant_fxn=MutantCategory.rmv_score_mutant,
            get_mutant_fxn=MutantCategory.rmv_get_mutant,
            add_muts_fxn=MutantCategory.meta_add_muts,
            motif_type=meta_motifs,
            max_iter=20,
            nts_per_iter=1,
            muts_to_keep=2,
            mut_iter_branch=10,
            num_feat_to_mut=1,
            min_msets_per_iter=10,
            score_key_fxn=lambda mc, new_score, orig_score: abs(new_score),
            call_fxn=MutantCategory.rmv_call_fxn,
            min_separation=3)

#===============================================================================
# conserved - MAMCONSERV FEATURES
#===============================================================================
#conserved 1 nt mutations
mut_cats['cnsrv_1nt'] = \
    MutantCategory(\
        cat='cnsrv_1nt',
        score_mutant_fxn=lambda mc, mut: 0,
        get_mutant_fxn=MutantCategory.conserved_get_mutant,
        pick_nt_fxn=MutantCategory.conserved_mut_loc,
        motif_type='MamConserv',
        max_iter=1,
        nts_per_iter=1,
        muts_to_keep=1,
        mut_iter_branch=0,
        num_feat_to_mut=3)

#conserved 3 nt mutations
mut_cats['cnsrv_3nt'] = \
    MutantCategory(\
        cat='cnsrv_3nt',
        score_mutant_fxn=lambda mc, mut: 0,
        get_mutant_fxn=MutantCategory.conserved_get_mutant,
        pick_nt_fxn=MutantCategory.conserved_mut_loc,
        motif_type='MamConserv',
        max_iter=1,
        nts_per_iter=3,
        muts_to_keep=1,
        mut_iter_branch=0,
        num_feat_to_mut=3)
#===============================================================================
# controls - CONTROL FEATURES
#===============================================================================

#random region
for rnd_region in ('exon', 'intron'):
    #rnd num is a tuple; first number is nt_per_iter, second is muts to keep
    for rnd_num in ((1, 2), (2, 2), (3, 1), (5, 1)):
        cat_name = 'rnd_{}_{}nt'.format(rnd_region, rnd_num[0])
        mut_cats[cat_name] = \
            MutantCategory(\
                cat=cat_name,
                score_mutant_fxn=lambda mc, mut: 0,
                get_mutant_fxn=lambda mc, sr: sr.get_features(mc.rnd_reg),
                pick_nt_fxn=MutantCategory.random_mut_loc,
                motif_type=None,
                max_iter=1,
                min_msets_per_iter=0,
                nts_per_iter=rnd_num[0],
                muts_to_keep=rnd_num[1],
                mut_iter_branch=0,
                num_feat_to_mut=2,
                rnd_reg=rnd_region,
                min_separation=rnd_num[0])

#===============================================================================
# aggr_* - AGGRESSIVE RANDOM MUTATION OF INTRONS/EXONS
#===============================================================================
for mtype in ('exon', 'intron', 'both'):
    cat_name = 'aggr_{}'.format(mtype)
    mut_cats[cat_name] = \
        MutantCategory(\
            cat=cat_name,
            get_mutant_fxn=lambda mc, sr: mc.aggr_get_mutants(sr),
            pick_nt_fxn=lambda self, sr, feat, count: \
                MutantCategory.aggr_mut_loc(self, sr, feat, count),
            motif_type=mtype,
            max_iter=1,
            nts_per_iter=1,
            muts_to_keep=2,
            mut_iter_branch=0,
            num_feat_to_mut=2,
            min_separation=1)

#===============================================================================
# *_splice - SPLICE SIGNAL FEATURES
#===============================================================================

spl_dict = {'d': 'me_splice_donor', #donor
             'a': 'me_splice_acceptor'}                            #acceptor

#for no splice, same splice, and consen splice, etc, give a tuple:
#    [0] function that generates splice mutants (pick_nt_fxn)
#    [1] function that tells whether or not to keep a mutant (keep_mut_fxn)
#    [2] how many mutants to make


class_dicts = {
    'no_spl': {
        'pick_nt_fxn': lambda self, sr, feat, count: \
            MutantCategory.splice_mut_loc(self, sr, feat, count, invert=True),
        'keep_mut_fxn': lambda mc, new_score, orig_score: True,
        'muts_to_keep': 1,
        'score_key_fxn': lambda mc, new_score, orig_score: new_score,
        },

    'same_splice': {
        'pick_nt_fxn': MutantCategory.splice_mut_loc,
        'keep_mut_fxn': lambda mc, new_score, orig_score: \
                               abs(orig_score - new_score) < 0.5,
        'muts_to_keep': 2,
        'score_key_fxn': lambda mc, new_score, orig_score: \
                                abs(orig_score - new_score)
        },

    'weak_spl': {
        'pick_nt_fxn': MutantCategory.splice_mut_loc,
        'keep_mut_fxn': lambda mc, new_score, orig_score: \
                               0 < new_score < 3 and \
                               abs(orig_score - new_score) > 0.5,
        'muts_to_keep': 2,
        'score_key_fxn': lambda mc, new_score, orig_score: \
                                util.range_dist(new_score, (0, 3)),
        'skip_mut_fxn': lambda mc, f, sr: f.qualifiers['evidence'] < 3
        },

    'strong_spl': {
        'pick_nt_fxn': MutantCategory.splice_mut_loc,
        'keep_mut_fxn': lambda mc, new_score, orig_score: \
                                   new_score > 10 and \
                                   abs(orig_score - new_score) > 0.5,
        'muts_to_keep': 2,
        'score_key_fxn': lambda mc, new_score, orig_score:-new_score,
        'skip_mut_fxn': lambda mc, f, sr: f.qualifiers['evidence'] > 10
        },
    }

#create category mutant objects from class dicts
for spl_site, mtype in spl_dict.items():
    for spl_class, sc_attribs in class_dicts.items():
        #generate category name
        cat_name = '{}_{}'.format(spl_class, spl_site)
        #create mutant category objects
        mut_cats[cat_name] = \
            MutantCategory(\
                cat=cat_name,
                keep_mutant_fxn=sc_attribs['keep_mut_fxn'],
                score_mutant_fxn=lambda mc, mut: \
                    motif.motif_types[mc.motif_type].score(mut.str)['scores'][0],
                get_mutant_fxn=lambda mc, sr: \
                    sr.get_features(mc.motif_type, correct=True),
                pick_nt_fxn=sc_attribs['pick_nt_fxn'],
                motif_type=mtype,
                max_iter=5,
                nts_per_iter=2,
                muts_to_keep=sc_attribs['muts_to_keep'],
                mut_iter_branch=4,
                num_feat_to_mut=2,
                score_key_fxn=sc_attribs['score_key_fxn'],
                skip_mut_fxn=sc_attribs['skip_mut_fxn'] \
                    if 'skip_mut_fxn' in sc_attribs
                    else lambda mc, f, sr: False)

#separate set of mutations to remove spurious splice signals
for spl_site, mtype in spl_dict.items():
        #generate category name
        cat_name = 'csplice_{}'.format(spl_site)
        #create mutant category objects
        mut_cats[cat_name] = \
            MutantCategory(\
                cat=cat_name,
                keep_mutant_fxn=lambda mc, new_score, orig_score: \
                    0 > new_score,
                score_mutant_fxn=lambda mc, mut: \
                    motif.motif_types[mtype].score(mut.str)['scores'][0],
                get_mutant_fxn=MutantCategory.csplice_get_mutant,
                pick_nt_fxn=lambda self, sr, feat, count: \
                    MutantCategory.splice_mut_loc(\
                        self, sr, feat, count, full=True),
                motif_type=mtype,
                max_iter=5,
                nts_per_iter=2,
                muts_to_keep=1,
                mut_iter_branch=4,
                num_feat_to_mut=2,
                score_key_fxn=lambda mc, new_score, orig_score: new_score)

#===============================================================================
# RBPmats - RNA binding protein motifs
#===============================================================================
mut_cats['RBPmats'] = \
    MutantCategory(\
        cat='RBPmats',
        keep_mutant_fxn=lambda mc, new_score, orig_score: \
            new_score == 0,
        score_mutant_fxn=MutantCategory.rbp_score_mutant,
        get_mutant_fxn=lambda mc, sr: mc.meta_get_mutant(sr,
                                                          mtype='RBPmats',
                                                          meta=False),
        motif_type='RBPmats',
        max_iter=3,
        nts_per_iter=3,
        muts_to_keep=1,
        mut_iter_branch=2,
        num_feat_to_mut=5,
        min_msets_per_iter=3)

#------------------------------------------------------------------------------
#MULTI_MUTANT CATEGORIES - they are put in the multi_mut_cats dict
#                          called with mutate_and_iterate()
#                          because they combine already made mutants, they are
#                          separated from mut_cats and are called after mut_cats

#===============================================================================
# p_aggr_intr - combines both above intronic randomizations
#===============================================================================

multi_mut_cats['p_aggr_intr'] = \
    MutantCategory(\
        cat='p_aggr_intr',
        motif_type='aggr_intron',
        max_iter=1,
        muts_to_keep=2,
        num_feat_to_mut=1,
        get_mutant_fxn=lambda mc, sr: [(0, 170)],
        pick_nt_fxn=MutantCategory.paired_mut_loc,
        nts_per_iter=1
    )
#===============================================================================
# p_aggr_intr - combines both above intronic randomizations
#===============================================================================

multi_mut_cats['p_weak_spl'] = \
    MutantCategory(\
        cat='p_weak_spl',
        motif_type=('weak_spl_a', 'weak_spl_d'),
        max_iter=1,
        muts_to_keep=2,
        num_feat_to_mut=1,
        get_mutant_fxn=lambda mc, sr: [(0, 170)],
        pick_nt_fxn=MutantCategory.paired_mut_loc,
        nts_per_iter=1,
        skip_mut_fxn=lambda mc, f, sr: (
            mut_cats['weak_spl_a'].skip_mut_fxn(
                sr.get_features('me_splice_acceptor',
                                correct=True).next(),
                sr
            ) or
            mut_cats['weak_spl_d'].skip_mut_fxn(
                sr.get_features('me_splice_donor',
                                correct=True).next(),
                sr
            )
        )
    )

multi_mut_cats['p_strong_spl'] = \
    MutantCategory(\
        cat='p_strong_spl',
        motif_type=('strong_spl_a', 'strong_spl_d'),
        max_iter=1,
        muts_to_keep=2,
        num_feat_to_mut=1,
        get_mutant_fxn=lambda mc, sr: [(0, 170)],
        pick_nt_fxn=MutantCategory.paired_mut_loc,
        nts_per_iter=1,
        skip_mut_fxn=lambda mc, f, sr: (
            mut_cats['strong_spl_a'].skip_mut_fxn(
                sr.get_features('me_splice_acceptor',
                                correct=True).next(),
                sr
            ) or
            mut_cats['strong_spl_d'].skip_mut_fxn(
                sr.get_features('me_splice_donor',
                                correct=True).next(),
                sr
            )
        )
    )
##-------------------------------------------------------------------------------
##OTHER MUTANT CLASSES - they go in the mutate.other_mutants dict and are not
##                       called via mutate_and_iterate, though their mutants
##                       can still be used by multi_mut_cats
#
##===============================================================================
## rmv_* - removes all instances of a motif, except in splice signals
##===============================================================================
#
#rmv_motif_cats = {}
#for mtype in motif.motif_types.keys():
#
#    multi_mut_cats['rmv_'+mtype] = \
#    MutantCategory(
#        cat = 'rmv_'+mtype,
#        keep_mutant_fxn = lambda mc, new_score, orig_score: \
#            abs(orig_score) * 0.1 > abs(new_score),
#        score_mutant_fxn = MutantCategory.meta_score_mutant,
#        get_mutant_fxn = lambda mc,sr: [],
#        motif_type = mtype,
#        max_iter = 5,
#        nts_per_iter = 2,
#        muts_to_keep = 1,
#        mut_iter_branch = 5,
#        num_feat_to_mut = 1,
#        min_msets_per_iter = 2,
#        call_fxn = lambda self,sr: remove_all_motifs(sr, self, self.motif_type))
#
##-------------------------------------------------------------------------------



def codon_fwd_table():
    '''generate a synonymous fwd-table for standard codons'''
    return CodonTable.standard_dna_table.forward_table

def codon_back_table():
    '''generate a synonymous back-table for standard codons'''

    #get Bio.Data's standard forward codon table
    fwdt = codon_fwd_table()

    #generate an empty backward table dict (codon -> set(AAs))
    bckt = {}
    for aa in set(fwdt.values()):
        bckt[aa] = set()

    for codontuple in fwdt.items():
        codon = codontuple[0]
        aa = codontuple[1]
        bckt[aa].add(codon)

    for aa in bckt.keys():
        bckt[aa] = frozenset(bckt[aa])

    return bckt

def mutant_locations(sr):
    #UPDATE TO BELOW: 10/24/11: I am skipping this and allowing two base
    #differences since I am checking for synonymous codons when I make mutants
    #to features....
    #
    #NOTE: One caveat of this method is that it does generate all one-base
    #mutants possible, but you can never reach synonymous codons that require
    #two base mutations, even if you combine multiple mutant_locations members
    #from the resulting muts array. I'm allowing this for now, since it will be
    #a  small percent of the mutation space.
    nucs = set('GATC')
    bckt = codon_back_table()
    fwdt = codon_fwd_table()

    exons = \
        interval(*[f.extract_pos() for f in sr.get_features('exon')])

    #separate exonic and intronic portions
    introns = interval((0, len(sr.seq))) & exons.invert()

    muts = {}

    #make non-coding mutant sets
    for intron_i in introns.components:
        for loc in range(*intron_i.to_tuple()):
            nuc = str(sr.seq[loc]).upper()
            muts[(loc, loc + 1)] = (nucs - set(nuc))

    #make coding mutant sets
    exp_exon_tup = expand_motif_to_codons(sr, exons.to_tuple()) #not necessary?
    codon_intervals = range(exp_exon_tup[0], exp_exon_tup[1], 3)

    for c_loc in codon_intervals:
        codon = str(sr.seq[c_loc:(c_loc + 3)]).upper()
        other_codons = bckt[fwdt[codon]]
        #see UPDATE above -> #pick codons that differ by exactly one base
        for other_cod in other_codons:
            diffs = util.str_diff(other_cod, codon)
            #if len(diffs) == 0 or len(diffs) > 1:
            #    continue
            for diff in diffs:
                diff_loc = (c_loc + diff, c_loc + diff + 1)
                if diff_loc in muts:
                    muts[diff_loc].add(other_cod[diff])
                else:
                    muts[diff_loc] = set(other_cod[diff])

    return muts

def add_muts_to_set(sr, mut, mut_cat,
                    nts_to_add=None,
                    num_new_sets=None,
                    seen_muts=frozenset(),
                    loc=None):
    '''
    need to update this docstring
    '''

    if num_new_sets == None:
        num_new_sets = mut_cat.new_per_branch
    if nts_to_add == None:
        nts_to_add = mut_cat.nts_per_iter

    seen_sets = frozenset([seen_mut.set for seen_mut in seen_muts])

    #remove positions that are already mutated from this set
    disj_ivl = interval(mut.loc) if loc == None else interval(loc)
    disj_ivl = disj_ivl & interval(*[tup[0] for tup in mut.set]).invert()

    #maximum new mutant sets to make
    max = mut_cat

    #this iterator will generate sets of new mutants that do not appear in
    # the current set
    feature.random_seed = hash(seen_sets) ^ hash(mut) ^ hash(disj_ivl)
    new_mut_addition_iter = sr.mutant_locations(disj_ivl, nts_to_add)

    #this counter will hold the number of mutant sets generated so far
    new_mut_sets_genned = 0

    #this list will hold the new mutant sets
    new_mut_sets = []

    #use the iterator to make new mutants and add them to the above list if they
    # are new
    for new_muts in new_mut_addition_iter:
        if new_mut_sets_genned >= num_new_sets: break

        new_mut = mut.with_new_mutants(new_muts)
        if new_mut not in seen_sets and new_mut.is_synonymous():
            new_mut_sets.append(new_mut)
            new_mut_sets_genned += 1

    #if there are no new mutants with $num muts added, then try with $num-1
    # recursively
    if (new_mut_sets_genned < num_new_sets
        and nts_to_add > 1):
        new_mut_sets.extend(
            add_muts_to_set(sr, mut, mut_cat,
                            nts_to_add=nts_to_add - 1,
                            num_new_sets=num_new_sets - new_mut_sets_genned))

    return new_mut_sets

def mutate_and_iterate(sr, feat, mut_cat, seed=12345, orig_score=None):
    '''
    psdueocode:

    1 get relevant variables and create object to sort mutant sets
    2 iterate forever (until max iter, criteria met, or dead end iter):
        2.1 convert mutant sets to strings
        2.2 for each set/string pair:
            2.2.1 score set/string with external fxn
            2.2.2 add score/set to list of pairs sorted by score
            2.2.3 trim the best list down to the max required
            2.2.3 after each set/str, if the entire best list matches criteria
                  then break
            2.2.4 also keep a branch list, potentially longer than best list
                  to carry forward to the next iteration, and trim it as well
        2.3 add more mutants to the branch sets with add_muts_to_set
            takes (sr, loc, mset, nts_per_iter, max_mut_sets, seen_sets)
        2.4 add old branch sets to this new list
        2.5 if this new list is same as one from prev iteration, then stop
            and give 'DeadEndIter'
    3 With our best sets (criteria met or not), return mutant objects

    TODO:
         - make 3 a separate function
         - make the loc passed in 2.3 be modifiable by an external fxn

    TO MAKE IT ITERABLE AT THE LEVEL OF MULTIPLE FEATURES
        - score function should take the sum of all features
        - loc passed to add_muts_to_sets should add muts only to features
        - 2.3 should be switchable to another function
          (namely, another mutate_and_iterate function)

    '''

    #parse inputs, seqfeature object versus tuple
    if isinstance(feat, tuple):
        loc = feat
    elif isinstance(feat, SeqFeature):
        loc = feat.extract_pos()
        try: orig_score = feat.extract_score()
        except: pass #take orig_score from kwargs

    #--------------------------------------------------
    #define relevant local variables for use during iteration
    #and initialize starting mutants
    max_iter = mut_cat.max_iter
    nts_per_iter = mut_cat.nts_per_iter
    min_msets_per_iter = mut_cat.min_msets_per_iter
    min_sep = mut_cat.min_separation
    mut_iter = (i for i in mut_cat.pick_nt_fxn(sr, feat, nts_per_iter))
    num_to_keep = mut_cat.muts_to_keep

    #define a function to extract a comparable set of built-in types
    get_froz_set = lambda mset: frozenset([frozenset(mset) for m in mset])

    mut_sets = []
    count = nts_per_iter
    for i in range(mut_cat.new_per_branch):
        try:
            mut_sets.append(mut_iter.next())
        except StopIteration:
            if count == 1 or len(mut_sets) >= min_msets_per_iter:
                break
            else:
                count -= 1
                mut_sets.extend(mut_cat.pick_nt_fxn(sr, feat, count))
    iter = 1

    #convert mutant sets into mutant objects (which also scores them)
    last_frozen_set = frozenset()
    muts = sortedset([Mutant(sr, mset,
                             loc=loc,
                             orig_score=orig_score,
                             category=mut_cat,
                             iter=iter) for mset in mut_sets],
                             key=lambda mut: mut.__key__())

    #keep track of previous version of the set
    #notes dict will be added to all muts
    notes = {}
    #--------------------------------------------------
    #enter our iterations of scoring & mutant_adding
    #break if we get the min_needed score or if we hit max_iter
    while True:
        #if all our SEPARATED mutants pass criteria, then we're good.
        #otherwise, keep looking
        #stop conditions----------------
        if iter >= max_iter:
            notes['MaxIterReached'] = True
            break
        if len(last_frozen_set) == len(muts):
            notes['NoNewMutants'] = True
            break
        if len(muts) == 0:
            notes['NoMutablePositions'] = True
            break
        sep_muts = list(Mutant.separated(muts, min_sep))
        if (len(sep_muts) >= num_to_keep
            and all(mut.criteria_met() for mut in sep_muts[:num_to_keep])):
            break
        #-------------------------------
        last_frozen_set = get_froz_set(muts)

        #MAKE BRANCH MUTS----------------------------------
        #to ensure that we sample a large enough space,
        #try to take mutants that are at least 2 mutations apart
        branchsize = min(len(muts), mut_cat.mut_iter_branch)
        branch_sets = Mutant.separated_best(muts, branchsize, min_sep)
        #ADD NEW NUCLEOTIDES-------------------------------
        #if we haven't found our mutant set, then add more mutants
        new_muts = \
            [mut_cat.add_muts_fxn(sr, m, seen_muts=muts) for m in branch_sets]
        muts.update(chain.from_iterable(new_muts))
        iter += 1

    #this returns muts regardless of separation (assuming we hit maxiter)
    sep_best = Mutant.separated_best(muts, num_to_keep, min_sep)
    [setattr(mut, 'feature', feat) for mut in sep_best
        if isinstance(feat, SeqFeature)]
    [mut.notes.update(notes) for mut in sep_best]
    return sep_best

def check_translation(dna1, dna2):
    '''check that two frame-0 dna strings have the same translation'''

    trans1 = Seq(dna1, generic_dna).translate().tostring()
    trans2 = Seq(dna2, generic_dna).translate().tostring()

    #print "{0} ({1}) vs\n{2} ({3})".format(trans1,dna1,trans2,dna2)
    if trans1 == trans2:
        return True
    else:
        return False


def expand_motif_to_codons(seq_record, motif_bounds):
    '''go through each exonic region and see if expanding the motif_bounds
    is necessary to include whole codons
    '''

    exon_list = seq_record.exon_list

    #generate a list of internal exon/intron start/ends within the motif
    (e_coords, i_coords) = get_motif_boundaries(motif_bounds, seq_record)

    #go through each exonic region and see if expanding the motif_bounds
    #is necessary to include whole codons
    for exonic_region in e_coords:
        exon_idx = find_enclosing_feature(exonic_region, exon_list)
        exon_bounds = exon_list[exon_idx].extract_pos()

        exp_exon = expand_region_to_codons(exon_bounds, exonic_region)

        if exp_exon[0] < motif_bounds[0]:
            motif_bounds = (exp_exon[0], motif_bounds[1])
        if exp_exon[1] > motif_bounds[1]:
            motif_bounds = (motif_bounds[0], exp_exon[1])

    return motif_bounds
#

def get_motif_boundaries(motif_slice, seq_record):
    '''
    return pythonic coordinate sets that are inside/outside of exon
     a tuple (exon/intron) of lists of tuples (coordinates within motif)
    '''

    #will put exonic/intronic overlapping coords into 2 lists
    exonic_overlap = []
    intronic_overlap = []

    exon_list = seq_record.exon_list


    #assume by default that sequence starts with exon, and cycle types as
    #exon/intron/exon/intron, but check this anyway and replace the cycle
    #if necessary

    #no exons exist
    if len(exon_list) == 0:
        exon_list = [FeatureLocation(0, len(seq_record.seq))]
        seq_type_cycle = cycle([intronic_overlap, exonic_overlap])

    #get list of exon starts/ends and chain them, so it looks like:
    #[exon1_start,exon1_end,exon2_start,exon2_end,...,seq_end]
    boundaries = list(chain(*map(lambda ex: ex.extract_pos(), exon_list)))

    #if the sequence does not end with an exon, append an intronic region to
    #the end
    #exons exist, but exon 1 does not start at 0
    if boundaries[0] != 0:
        seq_type_cycle = cycle([intronic_overlap, exonic_overlap])
        boundaries[:0] = [0]
    else:
        seq_type_cycle = cycle([exonic_overlap, intronic_overlap])


    if boundaries[-1:] != len(seq_record.seq):
        boundaries.append(len(seq_record.seq))

    #each index is a ref to the empty exonic/intronic coord lists
    seqtypes = list(
        islice(seq_type_cycle , len(boundaries) - 1))

    motif_srt = motif_slice[0]
    motif_end = motif_slice[1]

    #motif start is between boundaries[start_idx-1] and boundaries[start_idx]
    start_idx = bisect_left(boundaries, motif_srt)

    #motif end is between boundaries[end_idx-1] and boundaries[end_idx]
    end_idx = bisect_left(boundaries, motif_end)

    for region_idx in range(max(0, start_idx - 1), end_idx):
        try:
            reg_srt_idx = max(boundaries[region_idx], motif_srt)
            reg_end_idx = min(boundaries[region_idx + 1], motif_end)
        except IndexError:
            #if for some reason the bisection of the exon/intron boundaries with
            # the feature boundaries is not working...
            pass

        seqtypes[region_idx].append((reg_srt_idx, reg_end_idx))

    return (exonic_overlap, intronic_overlap)
#

def expand_region_to_codons(exon_coords, internal_coords):
    '''
    given an exonic coordinate set, expand the set to start at the first prior
    exon and end at the end of the next prior exon. This assumes that NO CODONS
    ARE SPLIT BETWEEN MULTIPLE EXONS.
    '''

    i_srt = internal_coords[0]
    i_end = internal_coords[1]
    e_srt = exon_coords[0]
    e_end = exon_coords[1]

    if i_srt < e_srt or i_end > e_end:
        raise(ValueError("Motif coordinates are not inside exon!"))

    #shift the i_srt back to start of first codon
    c_srt = i_srt - ((i_srt - e_srt) % 3)

    #shift the i_end to the end of last codon
    c_end = i_end + ((e_end - i_end) % 3)

    if c_end > e_end or c_srt < e_srt:
        raise(ValueError("This program can't deal with codons split between" +
                          "two exons, sorry!"))

    return (c_srt, c_end)

def find_enclosing_feature(coords, seq_feature_list):
    '''
    from an list of start/end tuples, find which one contains a given start/end
    tuple completely, returning the index.
    '''

    #negative one as placeholder index
    idx = -1

    #extract start and end position from each seq feature
    feat_coord_list = map(lambda sf: sf.extract_pos(), seq_feature_list)


    for feat_c_idx in range(len(feat_coord_list)):
        feat_c = feat_coord_list[feat_c_idx]
        if coords[0] >= feat_c[0] and coords[1] <= feat_c[1]:

            if idx != -1:
                warn('Multiple exons ' + str((idx, feat_c_idx)) +
                     ' contain coordinate set ' + str(coords) +
                     ', this is likely bad.')

            idx = feat_c_idx

    if idx == -1:
        raise(ValueError('Coordinates ' + str(coords) + ' were not found in' +
                         ' supplied coordinate list. Whoops.'))

    return idx

def overlaps_correct_motif(feat, seq_record):
    '''If this cryptic splice feature overlaps a good splice feature, we want
    to skip it, in case mutating it also affects real splice site. this will be
    dealt with more intelligently in the future.
    '''
    exon_list = seq_record.exon_list
    correct_ss = []

    for ss_feat in seq_record.features:
        if not feature.is_splicemod_motif(ss_feat):
            if 'function' not in ss_feat.qualifiers:
                continue

        if seq_record.is_correct_motif(ss_feat):
            correct_ss.append(ss_feat.extract_pos())



    retval = False
    fPos = feat.extract_pos()

    for cSS in correct_ss:
        if (cSS[0] <= fPos[0] <= cSS[1] or fPos[0] <= cSS[0] <= fPos[1]):
            retval = True

    return retval

def expand_mutants(mut_motifs, bounds, expand, seq_record):
    '''
    Given a list of mutant motif strings, it adds surrounding sequence and
    returns a new string list, move by expand[0] upstream and expand[1]
    downstream. Note that if you want to expand in both directions, expand[0]
    should be negative.
    '''

    seq_len = len(seq_record.seq)
    new_bounds = (bounds[0] + expand[0], bounds[1] + expand[1])

    #check to see if expansion will go outside of sequence bounds
    if new_bounds[0] < 0:
        expand = (0, new_bounds[1])
    if seq_len < new_bounds[1]:
         expand = (new_bounds[0], seq_len)

    prepend = seq_record.seq.tostring()[new_bounds[0]:bounds[0]]
    append = seq_record.seq.tostring()[bounds[1]:new_bounds[1]]
    return map(lambda seq: prepend + seq + append, mut_motifs)

def capitalize_mutant_string(orig_seq, mseq):
    ''' make the 'mutation' qualifier list nucleotides that have
        been changed in caps, like this: attacGacag
    '''
    orig_mut_letter_pairs = zip(list(orig_seq), list(mseq))
    return ''.join([m.lower() if o == m else m.upper() for o, m in
                    orig_mut_letter_pairs])

#def generate_feature_mutants(record,feat,new_scores,mask,
#                             values= cfg.gen_mut_count):
#    ''' arguments:
#        record:        a sequence record
#        feat:          the original feature to mutate
#        new_scores:    list of tuples of score ranges for mutants
#        mask:          a list of positions (ints) to keep constant when mutating
#        values (opt):  how many mutants to find for each new_score range in the
#                       list (defaults to 2)
#
#        This function will genereate len(new_scores) mutant features based on
#        the original feature. It will generate mutants with new scores based on
#        new_scores, which is a list of tuples containing the lowest and highest
#        allowable scores for each mutant. The mask consists of a list of
#        positions that will be invariant (such AGs).
#
#        It will update orig feat with a 'mutants' field that contains a dict
#        of dicts of tuples, a str and a float, that are sorted into keys that
#        correspond to new_scores.
#
#        eg: feat.qualifiers['mutants'] is a dict of dicts of tuples
#        eg: feat.qualifiers['mutants'][(3,4)] is a dict of tuples whose scores
#            are from 3 to 4
#        eg: feat.qualifiers['mutants'][(3,4)][(aagaTagat,3.0)] is a tuple that
#            corresponds to a mutant that differs by one nucleotide (T) and has
#            a score of 3.
#    '''
#
#    #test use: mutate.generate_feature_mutants(
#    #              record,feat,[(1.5,2.5),(3.5,4.5),(10,12)],[19,20])
#
#    feat.qualifiers['mutants'] = {}
#
#    #get feature location and sequence
#    feat_loc = feat.extract_pos()
#    orig_seq = record.seq[slice(*feat_loc)]
#    orig_score = feat.qualifiers['evidence']
#
#    def search_single_bp_mutants(record,feat,new_scores):
#        '''
#        This subfunction is to be called iteratively so we can search the
#        mutation space a single base at a time in order to find all mutants that
#        match our criteria.
#
#        Returns an updated version of the update_feat's mutant sequence
#        dictionary with the new mutants and their scores.
#        '''
#        #mutate sequences by one base pair
#        muta_seqs = mutate_motif(seq_feature= feat,seq_record= record, )
#
#        #score mutated sequences
#        muta_scores = score.call_max_ent( \
#                                feat.qualifiers['label'][0],muta_seqs)
#
#        #put each mutated seq into the qualifier dict in the original feature
#        #and check to ensure that each mutant is 'needed', i.e. it does not
#        #include a disallowed mutation and its score is in one of the ranges
#        #specified. If it meets these criteria, then we will add it and its
#        #score as a tuple to the feature's 'mutants' field in the qualifer dict.
#        for mseq, mscore in zip(muta_seqs,muta_scores):
#
#            needed_score = False
#            disallowed_mutation = False
#
#            #CHECK SCORE RANGES
#            inranges = [(lo,hi) if (mscore <= hi and mscore >= lo) else None \
#                        for lo,hi in new_scores]
#            inranges = filter(None,inranges)
#
#            if len(inranges) == 0: continue
#
#            mmutant = capitalize_mutant_string(orig_seq,mseq)
#
#            #CHECK MASK POSITIONS
#            #if a position index exists in mask, it is disallowed; skip it and
#            #don't include that mutant
#            orig_mut_letter_pairs = zip(list(orig_seq),list(mseq))
#            for disallowed_index in mask:
#                o,m = orig_mut_letter_pairs[disallowed_index]
#                if o != m: disallowed_mutation = True
#            if disallowed_mutation:
#                #print "SKIPPED: %s" % mmutant
#                continue
#
#            for inrange in inranges:
#                if inrange not in feat.qualifiers['mutants'] or \
#                    len(feat.qualifiers['mutants'][inrange]) < \
#                        cfg.gen_mut_stored_sbps:
#
#                    feat.qualifiers['mutants'].\
#                        setdefault(inrange,{})[mmutant] = mscore
#
#        return feat
#
#    #call the subfunction for the first time, looking one base away and update
#    #the original feature
#    feat = search_single_bp_mutants(record,feat,new_scores)
#
#    #use a little subfunction to check if all mutants have been found.
#    #we want to make sure that the dict has a list of all required ranges
#    #and that each range has at least $values mutants within it.
#    generate_feature_mutants.missing_muts = []
#
#    #when this function is called, check for missing mutants and if
#    #necessary, update the missing mutant list
#    def check_missing_mutants():
#        generate_feature_mutants.missing_muts = []
#        for inrange in new_scores:
#            if inrange not in feat.qualifiers['mutants'] or \
#                len(feat.qualifiers['mutants'][inrange]) < values:
#                generate_feature_mutants.missing_muts.append(inrange)
#        return len(generate_feature_mutants.missing_muts) > 0
#
#    #given a range, looks at all current mutants and returns a number of mutants,
#    #defaulting to 5, that are closest to the mutant score range, avoiding ones
#    #already tried
#    def get_mutants_near_to_range(range, old_muts,
#                                  num= cfg.gen_mut_follow_sbps): #MAGIC_NUMBER
#
#        #expand out all the range branches to get a list of tuples regardless of
#        #which range dict they are in
#        mut_tuples = list(
#            chain(*map(lambda r: r.items(),
#                       feat.qualifiers['mutants'].values())))
#
#        #filter out old muts
#        mut_tuples = [mt for mt in mut_tuples if mt[0].lower() not in old_muts]
#
#        r_siz = max(range)-min(range)/2
#        def get_range_dist(mut_tuple):
#            return max((0,util.range_dist(mut_tuple[1],range) - r_siz))
#
#        return sorted(mut_tuples, key=get_range_dist)[:(num)]
#
#    #if we haven't yet found at least 2 values for each score, iterate and keep
#    #trying progressively more mutated features until we have
#
#    #first, keep track of mutants that we've tried, as T/F dict of lowercase
#    #mut seqs
#    old_muts = {}
#
#    while check_missing_mutants():
#
#        #running check_missing_mutants will update the list of missing
#        #mutant ranges. We will now go through the missing ranges, find the
#        #two closest mutants to that score range, and mutate them to get
#        #closer to the range, hopefully hitting it.
#        close_mut_tuples = {}
#
#        #print "Missing mutants %s" % ' '.join(str(generate_feature_mutants.missing_muts))
#        #print "Qualifiers: %s" % str(feat.qualifiers)
#
#        for missing_range in generate_feature_mutants.missing_muts:
#
#            try_muts = get_mutants_near_to_range(missing_range,old_muts)
#            for near_mutants in try_muts:
#                close_mut_tuples[near_mutants] = True
#
#        #To mutate these mutants further, we have to mutate the record's raw
#        # sequence, but we will also change it back when we're done.
#        #It then adds these mutants to the original feat's mutant qualifier
#        # dict and cycles back to check for missing mutants using the while loop
#        for mutant_tuple in close_mut_tuples.keys():
#            mut_seq = mutant_tuple[0].lower()
#            #print "\tTrying %s..." % mutant_tuple[0]
#            record.replace_subseq(mut_seq,feat.extract_pos())
#            feat = search_single_bp_mutants(record,feat,new_scores)
#
#            #once we've searched for single bp derivatives starting from this
#            #sequence, save it in the array of old mutants so we don't do so
#            #again
#            old_muts[mut_seq] = True
#
#
#    record.replace_subseq(str(orig_seq),feat.extract_pos())
#
#    #print "Found all mutants: %s" % ' '.join(str(feat.qualifiers['mutants'].keys()))
#    #print "Qualifiers: %s" % str(feat.qualifiers)

def tups_to_str(seq_str, loc, tup_set):
    '''
    convert a set of tuples into a mutant string given a wild-type string
    '''
    mut_str = seq_str.upper()
    for tup in tup_set:
        #in the case of codon checks, the loc interval might not encompass
        #all of the mutant tuples. If so, then skip the mutants that are not
        #within the location given
        start = tup[0][0] - loc[0]
        end = tup[0][1] - loc[0]
        if start < 0 or end > len(mut_str): continue
        mut_str = ''.join([mut_str[:start], tup[1].lower(), mut_str[end:]])
    return mut_str

def str_to_tups(mut_str, orig_str, offset=0):
    '''
    given a mutant string and a wild type string, give a list of mutant tuples
    '''

    tups = ()

    diffs = util.str_diff(mut_str, orig_str, case_sensitive=False)

    for diff in diffs:
        diff_loc = (offset + diff, offset + diff + 1)
        diff_nt = mut_str[diff]
        tups += ((diff_loc, diff_nt),)

    return tups












