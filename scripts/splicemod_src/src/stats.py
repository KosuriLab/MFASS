'''
Created on Sep 15, 2011

@author: dbgoodman
'''
import numpy as np
from itertools import chain, ifilter
from collections import OrderedDict
import motif
import score
import util
import pdb


class SRStats:
    ''' This class defines a SeqRecord stats object, that holds statistical
        information for sequence records related to motifs, wiggle tracks,
        etc.
    '''

    #===========================================================================
    #GENERAL STATS
    #===========================================================================

    CALL_STATS = OrderedDict()

    #sr here is the seqrecord
    general = OrderedDict([
               ('chr',       lambda sr: sr.annotations['chr']),
               ('syn_range', lambda sr: '{}:{}'.format(\
                                         sr.annotations['synth_start'],
                                         sr.annotations['synth_end'])),
               ('id',        lambda sr: sr.id),
               ('up',        lambda sr: sr.annotations['synth_us']),
               ('exon_len',  lambda sr: sr.annotations['len']),
               ('down',      lambda sr: sr.annotations['synth_ds']),
               ('ss_num',    lambda sr: len(sr.correct_motifs)),
               ('ss_score',  lambda sr: score.score_seq_record(sr))
               ])

    CALL_STATS['general'] = ([''],
                             general,
                             lambda sr,x: (sr,))

    #===========================================================================
    #MOTIF STATS
    #===========================================================================
    #x here is a feature list iterator, sr is seq_record
    motif_attribs = OrderedDict([
        ('sum'  , lambda sr, x: sum((f.extract_score() for f in \
                                    sr.get_features(x.type)))),
        ('count', lambda sr, x: sum((1 for f in sr.get_features(x.type)))),
        ('comp' , lambda sr, x: sum((1 for f in sr.get_features(x.type,\
                                                               meta=True)))),
        ('ratio', lambda sr, x: sum((1 for f in sr.get_features(x.type))) / \
                               (sum((1 for f in sr.get_features(x.type,\
                                                                meta=True)))\
                                or float('-inf'))),
        ('max'  , lambda sr, x: util.max_none((f.extract_score() \
                                              for f in \
                                              sr.get_features(x.type,\
                                                              meta=True))))
        ])

    CALL_STATS['motif'] = (sorted(motif.motif_types.values()), # x values
                           motif_attribs, # stats to get (dict above)
                           lambda sr, x: (sr,x)) #how to get


    #===========================================================================
    #AGGREGATE STATS
    #===========================================================================
    #ag is aggregate group, sr is seq record
    agg_attribs = OrderedDict({
        'sum'  : lambda sr, ag: sum((f.extract_score() for f in \
                                     sr.get_features(ag,aggregate=True)))
        })

    CALL_STATS['agg'] = (motif.agg_types, # x values
                         agg_attribs, # stats to get (dict above)
                         lambda sr, ag: (sr, ag)) #how to get

    #===========================================================================
    #WIG STATS
    #===========================================================================
    #x here is a list of wig arrays for features, sr is seq features
    #TODO: MAKE THIS FOR ALL POTENTIAL TRACKS, NOT JUST MAMCONSERV
    wig_attribs = OrderedDict([
        ('cnsv_avg'   , lambda sr, x: util.mean_none([fl for fl in chain(*x)])),
        ('cnsv_f_max' , lambda sr, x: util.max_none(util.mean_none(f) for f in x))
        ])

    CALL_STATS['wig'] = (list(chain(sorted(motif.motif_types.values()),\
                               ('exon','intron'))),
                         wig_attribs,
                         lambda sr, x: (sr,sr.wigs_for_ftype(x,'MamConserv')))

    #===========================================================================
    #SNP STATS
    #===========================================================================

    #x here is a list of snps (variation features), sr is seq_record
    snp_attribs = OrderedDict([
        ('count' , lambda sr,x: sum((1 for f in x))),
        ('syn'   , lambda sr,x: sum(('synonymous_codon' in snp.qualifiers['note'] \
                                 for snp in x))),
        ('nonsyn', lambda sr,x: sum((has_effect(nonsyn_var,snp) \
                                 for snp in x))),
        ('intron', lambda sr,x: sum(('intron_variant' in snp.qualifiers['note'] \
                                 for snp in x))),
        ('splice', lambda sr,x: sum((has_effect(splice_var,snp) \
                                 for snp in x))),
        ('bad',    lambda sr,x: sum((has_effect(bad_effects,snp) \
                                 for snp in x)))
    ])

    CALL_STATS['snp'] = (['snp'], # x values (here only one)
                         snp_attribs, # stats to get (dict above)
                         lambda sr, x: (sr,sr.get_features('variation'),))


    def __init__(self,sr):

        self._stats_odict = OrderedDict()

        sr.stats = self

        for stat_group in SRStats.CALL_STATS:
            stat_actions = SRStats.CALL_STATS[stat_group]
            oloop_list, fxn_odict, get_st_fxn = stat_actions

            self._stats_odict[stat_group] = OrderedDict()

            for i in oloop_list:
                self._stats_odict[stat_group][str(i)] = OrderedDict()

                for st in fxn_odict:
                    lambda_args = get_st_fxn(sr,i)
                    self._stats_odict[stat_group][str(i)][str(st)] = \
                        fxn_odict[st](*lambda_args)

    def __getitem__(self,*args):
        return self._stats_odict.__getitem__(*args)

    def __setitem__(self,*args):
        return self._stats_odict.__setitem__(*args)

    def __iter__(self):
        return self._stats_odict.__iter__()

    def __str__(self):

        out_str = ''

        for sc in self:
            for st in self[sc]:
                for stat in self[sc][st]:
                    out_str += "{}\t".format((self[sc][st][stat]))

        return out_str

    def header(self):

        out_str = ''

        for sc in self:
            for st in self[sc]:
                for stat in self[sc][st]:
                    out_str += "_".join([sc,st,stat])+"\t"

        return out_str


    def __repr__(self):

        out_str = ''

        for sc in self:
            out_str += "{}=============\n".format(sc)
            for st in self[sc]:
                out_str += "\t{}:\n".format(st)
                for stat in self[sc][st]:
                    out_str += "\t\t{}:\t{}\n".format(stat,self[sc][st][stat])

        return out_str

#http://www.ensembl.org/info/docs/variation/index.html#consequences
bad_effects = ['complex_change_in_transcript',
               'stop_gained',
               'frameshift_variant']

splice_var =  ['splice_region_variant',
               'splice_acceptor_variant',
               'splice_donor_variant']

nonsyn_var =  ['inframe_codon_gain',
               'inframe_codon_loss',
               'non_synonymous_codon']

has_effect = lambda effects,snp: (bool(sum([ef in snp.qualifiers['note'] \
                                            for ef in effects])))
