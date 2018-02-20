from bx import wiggle
from Bio import Motif
from Bio.Seq import Seq
from Bio.Motif.Thresholds import ScoreDistribution
from collections import defaultdict
import re
import subprocess
import acora
import gzip
import glob
import itertools
import pdb
import pickle
import os

from numpy import array, where, ones

import cfg
import util
import sys

motif_types = {}
wiggle_tracks = {}

def me_splice_note_str(motif_bounds, string):
    boundary = abs(motif_bounds[0])

    return '|'.join([string[0:boundary].upper(),
                     string[boundary:].lower()])

def default_note_str(motif_bounds, string):
    return string.upper()

class SeqMotifType:
    '''
    This class stores information and pointers to function required to score
    and mutate certain motif types, like 3' and 5' splice sites, branch points,
    splice signals, etc, etc.
    '''

    def __hash__(self):
        return hash(self.type)

    def __eq__(self, other):
        return isinstance(other, SeqMotifType) and self.type == other.type

    def __init__(self,
                 type,
                 score_type,
                 upstr=None,
                 invariant=[],
                 score_dict={},
                 bounds=None,
                 filter_score=None,
                 note_str_func=default_note_str,
                 # mutant_check= _mutant_check
                 ** attribs
                 ):
        self.type = type
        self.upstr = upstr
        self.invariant = invariant
        self.filter_score = filter_score
        self.note_str = note_str_func
        self.bounds = bounds
        self.score_type = score_type
        # self.mutant_check = mutant_check

        #------------------------------------------------
        # FIRST: figure out what score type this motif is.
        #------------------------------------------------

        # aho-corasick search tree
        if self.score_type == 'acora':
            # we need a score dict from a file and a separate acora object
            if not isinstance(score_dict, str):
                raise NeedScoreDictFileException
            self.score_dict = \
                dict([record.split() for record in open(score_dict)])

            self.acora_tree = acora.AcoraBuilder(self.score_dict.keys()).build()
            self.score = self.acora

        # ternary search tree
        elif score_type == 'tst':
            # make our file of nmers and scores a tst object
            self.score_dict = tst.TST()
            tstmap = lambda tuple: self.score_dict.put(*tuple)
            map(tstmap, ([record.split() for record in open(score_dict)]))

            self.score = self.tst

        # max ent nmer score
        elif score_type == 'max_ent':
            self.score_dict = {}
            self.score = self.max_ent

            # open up an 'interactive' pipe to the maxent software for this motif
            programs = {'me_splice_donor':'score5', 'me_splice_acceptor':'score3'}

            self.command = cfg.programTemplate.substitute(path=cfg.maxEntPath,
                                        program=programs[self.type])

        # positon frequency matrix list
        elif score_type == 'pfm':
            print "Loading position frequency matrices from {}...".format(
                    self.type)
            pfm_glob = glob.glob(score_dict)
            name_pattern = re.compile('.*/(\w+).pfm')

            self.score_dict = {}

            if cfg.pickle_pfms and os.path.exists(cfg.pickle_pfms):
                motif_objs_dict = pickle.load(open(cfg.pickle_pfms, 'rb'))
                for motif_name, motif_obj in motif_objs_dict.items():
                    self.score_dict[motif_name] = motif_obj

            else:
                motif_objs_dict = {}

                for motif_file in pfm_glob:
                    motif_obj = Motif.read(open(motif_file), 'jaspar-pfm')
                    motif_name = name_pattern.match(motif_file).group(1)
                    print "\t{}...".format(motif_name)
                    motif_obj.name = motif_name
                    if len(motif_obj) > 7:
                        self.score_dict[motif_name] = motif_obj
                        motif_obj.sd = \
                            ScoreDistribution(motif_obj, precision=10 ** 3)

                        # low false-positive rate to make sure motifs are real
                        motif_obj.thresh = max(1, motif_obj.sd.threshold_fpr(0.01))

                if cfg.pickle_pfms:
                    pickle.dump(motif_objs_dict, open(cfg.pickle_pfms, 'wb'))

            self.score = self.pfm
            print "Motif matrices done.\n"




        #------------------------------------------------
        # SECOND: parse filter score information.
        #------------------------------------------------

        # if filter score is an int, make it a lambda function that determines
        # whether or not it is a 'worthwile' score; this could be as simple as
        # > 0, or it could be a range, etc, etc. The lambda function will
        # return true if the score should be kept and false if it should not.

        if isinstance(self.filter_score, float) or \
            isinstance(self.filter_score, int):

            self.filter_score = \
                lambda val, min = self.filter_score: \
                    val > min

        elif isinstance(self.filter_score, tuple):

            self.filter_score = \
                lambda val, minmax = self.filter_score: \
                    val < minmax[0] or val > minmax[1]

        # if filter_score is none, always return true
        elif self.filter_score == None:
            self.filter_score = lambda val: True

        #------------------------------------------------
        # THIRD: add motif to motif type dict and cleanup
        #------------------------------------------------
        motif_types[self.type] = self
        self.attribs = attribs

    def __str__(self):
        return self.type

    @util.memoized
    def acora(self, string):
        '''
        this score function takes the already instantiated acora tree and
        returns matches
        '''

        nmers = []
        scores = []
        locations = []

        str_loc_tuples = self.acora_tree.findall(string.upper())
        for str, loc in str_loc_tuples:
            locations.append((loc, loc + len(str)))
            scores.append(self.score_dict[str])
            nmers.append(str)

        return {'nmers': nmers, 'locations': locations, 'scores': scores}

    @util.memoized
    def max_ent(self, string):
        '''
        score function for motifs that are found by maxent (splice acceptors
        and donors. like all score functions, returns nmers, locations, scores
        '''
        # print "Pre-score Subprocess state for {}: {}".format(self.type,self.p.poll())
        nmers, locations = util.generate_nmers(string, self.bounds)
        scores = self.call_max_ent(nmers)
        # print "Post-score Subprocess state for {}: {}".format(self.type, self.p.poll())
        # print "\tDid we get scores? {}".format(len(scores) > 0)


        return {'nmers': nmers, 'locations': locations, 'scores': scores}

    @util.memoized
    def pfm(self, string):

        nmers = []
        scores = []
        locations = []
        names = []

        for m_name, m in self.score_dict.items():
            for pos, score in m.search_pwm(Seq(string),
                                          threshold=m.thresh,
                                          both=False):
                locations.append((pos, (pos + len(m))))
                scores.append(score)
                names.append(m.name)
                nmers.append(string[pos:(pos + len(m))])

        return {'nmers': nmers, 'locations': locations, 'scores': scores}

    def call_max_ent(self, nmer_list):
        '''
        calls maxEnt programs given a list of n-mers. n-mers must be the right
        length
        '''
        nmer_list = array(nmer_list)

        scores = float('-inf') * ones(len(nmer_list))
        new_nmer_idxes = []

        # first, fill in precomputed scores, or mark new nmers for maxEnt run
        for nmer_i in range(len(nmer_list)):
            if nmer_list[nmer_i] in self.score_dict:
                scores[nmer_i] = self.score_dict[nmer_list[nmer_i]]
            else:
                new_nmer_idxes.append(nmer_i)

        if len(new_nmer_idxes) > 0:

            new_nmers_string = "\n".join(nmer_list[array(new_nmer_idxes)])

            self.p = subprocess.Popen(("perl", self.command + "_stdin.pl", "-"),
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   close_fds=False,
                                   cwd=cfg.maxEntPath,
                                   bufsize= -1)

            # line after this comment is broken, fixed below
            stdout_text, stderr_text = self.p.communicate(new_nmers_string)

            # split the output into an array
            new_scores = array(stdout_text.splitlines())
            # convert blank lines to -inf (they weren't long enough  to count)
            new_scores = where(new_scores == '', '-inf', new_scores)
            # convert to float
            new_scores = map(float, new_scores)
            # save new scores to old score list
            for new_nmer_i in range(len(new_nmer_idxes)):
                new_nmer = nmer_list[new_nmer_idxes[new_nmer_i]]
                if new_nmer not in self.score_dict:
                    self.score_dict[new_nmer] = \
                        new_scores[new_nmer_i]

            scores[new_nmer_idxes] = new_scores

        return scores

class Wiggle:
    '''
    Wiggle files (http://genome.ucsc.edu/goldenPath/help/wiggle.html) are files
    contain dense continuous genomic annotation data, like GC content, probabil-
    ity scores and transcriptome data. This class takes a set of wig files or
    or a gzipped set and loads them into iterators. Because it does not store
    then in memory, it is best that the function that uses a wiggle object does
    so sequentially throughout a chromosome instead of jumping around.

    Because this has been slow, I have added a seek feature using an index file
    for each wiggle file. I made these indices with a simple
    bash one-liner.

    The index file looks like:

    line_number \t byte_offset \t chr \t chr_pos_start

    where the first column is the position of a header line, such as:
    fixedStep chrom=chr1 start=34045 step=1

    the second column is the byte offset before said line starts, for mmap

    the third column is the chr position (in this case 34045) and the third
    line is the line number for this position.

    the script is:

    for wig in *.gz
    do
        gzcat $wig \
        | grep -nb fixedStep \
        | perl -ne 's/(\d+):(\d+).*start=(\d+).*/$1\t$2\t$3\t$4/ && print' \
        > $wig.idx
    done
    '''

    def __init__(self, name, wiggle_dir, file_prefix,
                 file_suffix, max_buffer=50):

        self.name = name

        self.lookback_buffers = {}
        self.max_buffer = max_buffer

        glob_components = (wiggle_dir, '/', file_prefix, '*', file_suffix)

        # put all files into a list
        gzlist = glob.glob(''.join(glob_components))

        # associate files into a dict by chromosome name
        re_components = map(util.to_raw, glob_components)
        # don't use wig files with '_'s because they are short contig
        # and there are too many
        re_components[-2] = "([a-zA-Z0-9]+)"

        pattern = re.compile(''.join(re_components))

        self.fnames = \
            dict([(pattern.match(gzfile).group(1), gzfile) \
                    for gzfile in gzlist if pattern.match(gzfile)])

        assert len(self.fnames) > 0, (
                'No wiggle files found in {}'.format(file_prefix))

        self.idxnames = \
            dict([(pattern.match(gzfile).group(1),
                    os.path.splitext(gzfile)[0] + r'.idx') \
                    for gzfile in gzlist if pattern.match(gzfile)])
        self.fhandles = \
            dict([(chr, gzip.open(fn)) \
                    for chr, fn in self.fnames.items()])
        self.fidxhandles = \
            dict([(chr, open(fn)) \
                    for chr, fn in self.idxnames.items()])

        # load the indices into a dict of chrs
        self.chridx = {}
        for chr, fidxh in self.fidxhandles.items():
            self.chridx[chr] = [map(int, fl.split()) for fl in fidxh.readlines()]

        wiggle_tracks[name] = self

    def get_region(self, chr, start, end):
        ''' we look through the chromosome index list and we pull out the
            seek position we want by finding the closest chr pos (col 3) that
            does not go over our start. we check to make sure that our end
            is less than the next chr pos (col 3 of next row). col 2 of the
            row is our seek() position. Then we figure out how many lines
            we need to skip, and how many we need to take, and we just return
            them all. no lookback buffer.

        chrom 0 , start 1, end 2, strand 3, value 4.
        '''

        start += 1  # to adjust for pythonic numbering
        end += 2

        next_idx_line = 0
        while self.chridx[chr][next_idx_line][2] < start:
            next_idx_line += 1

        if next_idx_line < 1: raise ValueError("START not in wiggle file")

        curr_idx_line = next_idx_line - 1

        next_idx_row = self.chridx[chr][next_idx_line]
        curr_idx_row = self.chridx[chr][curr_idx_line]

        # if the start and end positions are not on the same set of wiggles
        if next_idx_row[2] <= end:
            raise ValueError("This region is not contiguous in the wig file!")

        # if there are not enough lines in this wiggle set
        if next_idx_row[0] - curr_idx_row[0] < end - curr_idx_row[2]:
            import pdb
            pdb.set_trace()
            raise ValueError("This wiggle set is not large enough!")

        # get seek position, number of lines to jump, and number of lines to keep
        seekpos = curr_idx_row[1]
        jumplines = start - curr_idx_row[2]
        keeplines = end - start

        # tups to keep: chrom 0 , start 1, end 2, strand 3, value 4.
        tups = []

        fh = self.fhandles[chr]
        fh.seek(seekpos)
        for throwaway in range(jumplines):
            fh.readline()

        for keep in range(keeplines):
            # the start/end coords will be pythonic
            tups.append((chr, start + keep - 1, start + keep, '+', float(fh.readline())))

#        if chr not in self.lookback_buffers:
#            self.lookback_buffers[chr] = []
#
#        buffered_iter = \
#            itertools.chain(self.lookback_buffers[chr],self.wigiters[chr])
#        until_start = \
#            itertools.dropwhile(lambda x: x[1] < start - 1,buffered_iter)
#        until_end = \
#            itertools.takewhile(lambda x: x[1] < end, until_start)
#
#        map(tups.append, until_end)
#
#        if len(tups) != end - start + 1:
#            raise WiggleBadRangeSizeException
#
#        to_append = self.max_buffer - len(self.lookback_buffers[chr])
#
#        self.lookback_buffers[chr].extend(tups[-to_append:])

        return tups

#===============================================================================
# MOTIF TYPE DECLARATIONS
#===============================================================================

me_splice_donor = SeqMotifType(type='me_splice_donor',
                            score_type='max_ent',
                            upstr='intron',
                            invariant=[3, 4],
                            bounds=[-3, 6],
                            filter_score=1,
                            note_str_func=me_splice_note_str,
                            fiveprime='exon'
                            )

me_splice_acceptor = SeqMotifType(type='me_splice_acceptor',
                            score_type='max_ent',
                            upstr='intron',
                            invariant=[18, 19],
                            bounds=[-20, 3],
                            filter_score=1,
                            note_str_func=me_splice_note_str,
                            fiveprime='intron',
                            )

Ke2011_ESS = SeqMotifType(type='Ke2011_ESS',
                            upstr=None,
                            invariant=[],
                            score_type='acora',
                            score_dict=os.path.join(cfg.motifDir,
                                         'Ke2011/ESSseq.txt'),
                            filter_score=lambda val: float(val) < 0,
                            can_consolidate=True,
                            context='exon',
                            aggregate='Ke2011',
                            no_print=True
                            )

Ke2011_ESE = SeqMotifType(type='Ke2011_ESE',
                            upstr=None,
                            invariant=[],
                            score_type='acora',
                            score_dict=os.path.join(cfg.motifDir,
                                         'Ke2011/ESEseq.txt'),
                            filter_score=lambda val: float(val) > 0,
                            can_consolidate=True,
                            context='exon',
                            aggregate='Ke2011',
                            no_print=True
                            )
Vlkr07_DICS = SeqMotifType(type='Vlkr07_DICS',
                            upstr=None,
                            invariant=[],
                            score_type='acora',
                            score_dict=os.path.join(cfg.motifDir,
                                         'Voelker07/DI_CS.txt'),
                            can_consolidate=True,
                            context='donor_intron',
                            # aggregate= 'Ke2011',
                            no_print=True
                            )
Vlkr07_AICS = SeqMotifType(type='Vlkr07_AICS',
                            upstr=None,
                            invariant=[],
                            score_type='acora',
                            score_dict=os.path.join(cfg.motifDir,
                                         'Voelker07/AI_CS.txt'),
                            can_consolidate=True,
                            context='acceptor_intron',
                            # aggregate= 'Ke2011',
                            no_print=True
                            )
RBPmats = SeqMotifType(type='RBPmats',
                            upstr=None,
                            invariant=[],
                            score_type='pfm',
                            score_dict=os.path.join(cfg.motifDir,
                                         'jaspar/*.pfm'),
                            score_dist=os.path.join(cfg.motifDir,
                                         'jaspar/dists.txt')
                            )

# aggregate type explicit declarations
agg_types = ['Ke2011']

#===============================================================================
# WIGGLE DECLARATIONS
#===============================================================================

MamConserv = Wiggle(name='MamConserv',
                        wiggle_dir=os.path.join(
                                cfg.intronDataDir, 'hg19.100way.phyloP100way'),
                        file_prefix='chr',
                        file_suffix='.phyloP100way.wigFix.gz')








