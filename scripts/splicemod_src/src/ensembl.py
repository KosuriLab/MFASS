import MySQLdb
from interval import interval
from cogent.db.ensembl import HostAccount, Genome
# import cProfile

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

import ast
import random
import glob
import pdb

import motif
import cfg
import stats

ens_release = cfg.ens_release
db = MySQLdb.connect(**cfg.ens_rmt_db_dict)

# pycog = HostAccount(cfg.ens_lcl_db_dict['host'],
#                    cfg.ens_lcl_db_dict['user'],
#                    cfg.ens_lcl_db_dict['passwd'],
#                    cfg.ens_lcl_db_dict['port'])

pycog = HostAccount(cfg.ens_rmt_db_dict['host'],
                    cfg.ens_rmt_db_dict['user'],
                    cfg.ens_rmt_db_dict['passwd'],
                    cfg.ens_rmt_db_dict['port'])

hs37 = Genome('human', Release=78, account=pycog)
# hs37_remote = Genome('human', Release=63, account=pycog_remote)


def get_region(row_dict):
    return map(row_dict.get, ['chr', 'start', 'end'])

def type_me(str):
    try: return ast.literal_eval(str)
    except: return str

def categorize_exons(synth_size=cfg.chip_synth_length,
                     us_min=cfg.chip_us_min,
                     ds_min=cfg.chip_ds_min,
                     intron_padding=cfg.chip_intron_padding,
                     number_to_mutate=cfg.chip_number_to_mutate,
                     max_spots=cfg.chip_feature_count):
    '''
    for each gene:
    1. get all CCDS transcripts (if no CCDS transcripts, then skip)
    2. get all exons from those transcripts
    3. mark exons with a flag that meet following criteria:
       (a. between 12 and N bases long <-- skipping this for now)
       b. start and end on phase 0
    4. determine the size of region to synthesize, at least 40 bases up, 30 dwn
    5. based on this info, go through ALL exons (not just ccds), and for
       the following categories, fill in ccds_status, start, end:

       exons that overlap this exon
       exons that subsume this exon
       exons that this exon subsumes
       exons that have the same start and end
       exons that invade this exon's synthesized region
       exons whose flanks (50 bases) invade this exon's synthesized region

    '''

    # open files
    exon_file = open(cfg.ens_exon_fn)
    stats_file = open(cfg.ens_exon_stats, 'w')

    # set up counts dict
    counts = dict()
    counts['skip_size'] = 0
    counts['valid'] = 0
    counts['invalid'] = 0
    counts['total'] = 0
    counts['skip_identical_subsumed'] = 0
    counts['skip_overlap'] = 0
    counts['cut_sites'] = 0
    counts['dup_seq'] = 0

    # keep track of all 170mers

    seen_seqs = set()

    skip_exons = dict()  # list of exons to skip due to copies or being subsumed
    overlap_exons = dict()  # list of exons to skip due to being overlapped

    exon_dict_keys = exon_file.readline().split()

    all_exons = exon_file.readlines()
    # get random sample to perform mutations on
    random.seed(1)
    sample_indices = sorted(random.sample(range(len(all_exons)),
                            min(len(all_exons), number_to_mutate)))

    # \/ if we want to get them from mySQL
    #for i, ccds_ex in enumerate(get_ccds_exons(synth_size - us_min - ds_min)):


    # \/ if we want to get them from the ens_exon_fn file cache (in cfg.py)
    for i, ccds_ex in enumerate(all_exons):

        #------------------------------------------------------------------
        # 0. Keep time, running count, print every 50 exons
        #------------------------------------------------------------------

        counts['total'] += 1

        if counts['total'] > 0 and not (counts['total'] % 50):
            print "Finished %d exons..." % (counts['total'])
            for key in counts.keys():
                print "\t%s:\t%d" % (key, counts[key])


        #------------------------------------------------------------------
        # 1. Prepare CCDS_Ex exon record
        #------------------------------------------------------------------

        # if we are taking records from a file instead of mysql, then
        # split,map,zip into a dict
        ccds_ex = map(type_me, ccds_ex.split())
        ccds_ex = dict(zip(exon_dict_keys, ccds_ex))

        ccds_ex['invalid'] = 0
        ccds_ex['class'] = []
        ccds_ex['subsumed'] = []
        ccds_ex['identical'] = []

        #------------------------------------------------------------------
        # 2. Skip if this exon has been seen before and is
        #     subsumed, overlapped, or identical
        #------------------------------------------------------------------

        #if 'ENSE00001362516' not in ccds_ex['exon']: continue

        for fg in glob.glob(cfg.ens_fas_dir + '*'):
            if ccds_ex['exon'] in fg:
                print ccds_ex['exon'] + "\tPRESENT"
                fas_fn = (cfg.ens_fas_dir +
                         "{}/{}.fas".format(ccds_ex['exon'],
                                            ccds_ex['exon']))
                seq_record = SeqIO.parse(fas_fn, 'fasta').next()
                if str.upper(seq_record.seq.tostring()) not in seen_seqs:
                    seen_seqs.add(str.upper(seq_record.seq.tostring()))
                else:
                    print "\tDUPLICATE"
                skip_exons[ccds_ex['exon']] = 1


        if ccds_ex['exon'] in skip_exons:
            print ccds_ex['exon'] + "\tIDENTICAL"
            counts['skip_identical_subsumed'] += 1
            continue

        if ccds_ex['exon'] in overlap_exons:
            counts['skip_overlap'] += 1
            continue

        #------------------------------------------------------------------
        # 3. Get Region To Synthesize
        #------------------------------------------------------------------

        # synth_size is region to make, not counting primers/RE sites
        # us_min and ds_min are upstream and downstream minimums

        # skip if the exon plus us_min and ds_min is too large
        min_synth_size = ccds_ex['len'] + us_min + ds_min
        if min_synth_size > synth_size:
            counts['skip_size'] += 1
            continue

        synth_remainder, odd = divmod(synth_size - min_synth_size, 2)

        # if the exon is on strand -1, then we need to switch us and ds
        if ccds_ex['strand'] == 1:
            synth_us = us_min + synth_remainder + odd
            synth_ds = ds_min + synth_remainder
        else:  # opposite strand, switch us and ds
            synth_us = ds_min + synth_remainder
            synth_ds = us_min + synth_remainder + odd

        synth_region = dict()
        synth_region['chr'] = ccds_ex['chr']
        synth_region['start'] = ccds_ex['start'] - synth_us - 51
        synth_region['end'] = ccds_ex['end'] + synth_ds + 51

        ccds_ex['synth_start'] = ccds_ex['start'] - synth_us
        ccds_ex['synth_end'] = ccds_ex['end'] + synth_ds
        ccds_ex['synth_us'] = synth_us
        ccds_ex['synth_ds'] = synth_ds

        #------------------------------------------------------------------
        # 3. Find Exons that Exist in This Region
        #------------------------------------------------------------------

        for re_ex in get_exons_in_region(*get_region(synth_region)):

            re_ex['invalid'] = False
            re_ex['class'] = 'error'

            # region exon (this loop)
            re_ivl = interval[re_ex['start'], re_ex['end']]
            # ccds exon (outer loop exon)
            ce_ivl = interval[ccds_ex['start'], ccds_ex['end']]
            # synthesized region
            sr_ivl = interval[synth_region['start'], synth_region['end']]
            # flanking synthesized region
            fr_ivl = sr_ivl + interval(-50, 50)

            isiz = lambda ivl: ivl[0][1] - ivl[0][0]

            # look for situations that we want to avoid
            #3.1 CCDS Identical=================================================
            # exons that have the same start and end
            if re_ivl == ce_ivl:
                re_ex['class'] = 'identical'
                skip_exons[re_ex['exon']] = 1
                ccds_ex['identical'].append(re_ex['exon'])
                re_ex['invalid'] = False
            #3.2 CCDS Subsumes==================================================
            # exons that this exon subsumes
            # starts and ends inside of exon
            elif re_ivl in ce_ivl:
                re_ex['class'] = 'subsumes'
                skip_exons[re_ex['exon']] = 1
                ccds_ex['subsumed'].append(re_ex['exon'])
                re_ex['invalid'] = False
            #3.3 CCDS is Subsumed===============================================
            # exons that subsume this exon:
            #  starts and ends outside of exon
            elif ce_ivl in re_ivl:
                re_ex['class'] = 'subsumed'
                re_ex['invalid'] = True
            #3.4 Exon Overlap===================================================
            # exons that overlap this exon:
            # starts in middle, ends after or ends in middle, starts after
            elif (re_ivl[0][0] in ce_ivl and re_ivl[0][1] > ce_ivl[0][1])  \
               or (re_ivl[0][1] in ce_ivl and re_ivl[0][0] < ce_ivl[0][0]):
                re_ex['class'] = 'overlapped'
                re_ex['invalid'] = True
                overlap_exons[re_ex['exon']] = 1
            #3.5 Invasion=======================================================
            # exons that invade this exon's synthesized region
            elif re_ivl[0][0] in sr_ivl or re_ivl[0][1] in sr_ivl:
                re_ex['class'] = 'invaded'
                re_ex['invalid'] = True
            #3.6 Flank Invasion=================================================
            # exons whose flanks (50 bases) invade this exon's
            # synthesized region
            elif re_ivl[0][0] in fr_ivl or re_ivl[0][1] in fr_ivl:
                re_ex['class'] = 'flank_invaded'
                re_ex['invalid'] = True
            #===================================================================

            # FINALLY: outer loop exon invalid if subsumed, overlapped, invaded
            if re_ex['invalid']:
                ccds_ex['invalid'] = True
                ccds_ex['class'].append(re_ex['class'])

        # update class counts
        for exclass in ccds_ex['class']:
            if exclass in counts:
                counts[exclass] += 1
            else:
                counts[exclass] = 1

        # skip ccds exon if invalid
        if not ccds_ex['invalid']:
            counts['valid'] += 1
            print "%s\tOK" % ccds_ex['exon']
        else:
            counts['invalid'] += 1
            print "%s\t%s" % (ccds_ex['exon'], ccds_ex['class'])
            continue

        try:
            # only annotate & mutate if this exon is in the sample_indices set
            mutate_this = i in sample_indices

            #-------------------------------------------------------------------
            # 4. Get Seq Features and SNPs
            #-------------------------------------------------------------------
            # finally, make our exon into a bonafide splicemod SeqRecord
            exon_record = make_seq_record(ccds_ex,
                                          skip_annotation=not mutate_this)

            # check for cut sites
            for cs in cfg.cut_sites:
                if cs in exon_record.seq:
                    print " HAS {}".format(cs)
                    counts['cut_sites'] += 1
                    raise ValueError

            # check if an identical exon has been already seen
            if str.upper(exon_record.seq.tostring()) in seen_seqs:
                print " DUPLICATED"
                counts['dup_seq'] += 1
                raise ValueError

            #-------------------------------------------------------------------
            # 5. Make a list of mutants from MutantCategories
            #-------------------------------------------------------------------
            # make an attribute to store mutants in this seq record
            exon_record.mutants = {}
            if mutate_this:
                #5.1 Get Statistics=============================================
                exon_record.stats = stats.SRStats(exon_record)
                stats_str = str(exon_record.stats)
                # if this is the first stats str, print the header first
                if stats_file.tell() == 0:
                    stats_file.write(exon_record.stats.header() + "\n")
                stats_file.write(stats_str + "\n")
                #5.2 Generate Mutants===========================================
                exon_record.generate_all_mutants()

            #------------------------------------------------------------------
            # 6. Save to FASTA and GBK formats
            #------------------------------------------------------------------

            # this writes all of our mutants to genbank and fasta files
            exon_record.save_all_mutants(gb=mutate_this, fas=True)
            seen_seqs.add(str.upper(exon_record.seq.tostring()))

        except ValueError:
            print "EXON {} SKIPPED".format(exon_record.id)

    #---------------------------------------------------------------------------
    # Print final counts at the end:
    print "Final Counts:\n"
    for key in counts.keys():
        print "F\t%s:\t%d" % (key, counts[key])


def add_wiggle_data(exon_record):
    ''' sequentially add wiggle track information to exons. Because I am reading
        through chromosomes in order, I can just add wiggle data as an input
        stream using the wiggle class and the get_region method.
    '''
    chr = str(exon_record.annotations['chr'])
    start = exon_record.annotations['synth_start']
    end = exon_record.annotations['synth_end']

    for wig_tr in motif.wiggle_tracks.values():
        values = wig_tr.get_region(chr, start, end)
        if exon_record.annotations['strand'] == -1:
            values = list(reversed(values))

        if len(values) != cfg.chip_synth_length:
            raise ValueError("Wiggle track not long enough")

        exon_record.add_wiggle_track(wig_tr.name, map(lambda v: v[4], values))


def make_seq_record(exon, skip_annotation=False):
    ''' use pycogent to grab the sequence for the exon's region, and create a
        new seq record with all the bells and whistles. add in the predicted
        exon start and ends and splice sites as well.
    '''
    # GET SEQ FROM PYCOGENT
    #===========================================================================
    cog_reg = hs37.getRegion(CoordName=exon['chr'],
                                 Start=exon['synth_start'] - 1,
                                   End=exon['synth_end'])
    cog_seq = cog_reg.Seq

    # reverse complement if necessary
    if exon['strand'] == -1:
        cog_seq = cog_seq.rc()
    seq = cog_seq._seq

    if not skip_annotation:
    # DO SNPS
    #===========================================================================
        var = hs37.getFeatures(feature_types='variation', region=cog_reg)
        exon['snps'] = []
        for snp in var:
            if snp.Alleles != 'HGMD_MUTATION':
                exon['snps'].append(snp)

    # TODO: APPEND OUTER INTRON CONTEXT
    #===========================================================================
    # use the cfg information to get sequence, how much to add

    # MAKE SEQ RECORD
    #===========================================================================
    # create the seqRecord object, adding the exon dict as annotation fields
    record = SeqRecord(Seq(seq, generic_dna), id=exon['exon'],
                       annotations=exon)

    # add the exon features and putative splicing features, set the source as
    # ensembl_exon (not splicemod, so that it doesn't get deleted)

    if exon['strand'] == 1:
        in_front = exon['synth_us']
    else:
        in_front = exon['synth_ds']

    record.features = [\
        SeqFeature(FeatureLocation(in_front,
                   (exon['len'] + in_front)),
                   type="exon"),
        SeqFeature(FeatureLocation(0, in_front),
                   type="intron"),
        SeqFeature(FeatureLocation(in_front + exon['len'],
                   len(seq)),
                   type="intron")]

    for feat in record.features:
        feat.qualifiers['source'] = 'ensembl_exon'

    record.description = \
        ("{r.id} chr{chr}:{synth_start}-{synth_end} " +
         "strand={strand} len={synth_us}.{len}.{synth_ds} " +
         "ccds={CCDS}").format(r=record, **record.annotations)

    if not skip_annotation:
        record.populate_attribs()
        add_wiggle_data(record)
        record.add_conservation_features()

    return record

def get_exons_in_region(chr, region_start, region_end):
        exon_region_query = \
        '''
        select distinct
           ex.stable_id as 'exon',
           sr.name as 'chr',
           sr.seq_region_id,
           ex.seq_region_start as 'start',
           ex.seq_region_end as 'end'
        from
            exon as ex,
            seq_region as sr,
            transcript as tr,
            exon_transcript as extr
        where
            sr.name = '{0}' and
            sr.seq_region_id = ex.seq_region_id
            and ((ex.seq_region_start > {1} and ex.seq_region_start < {2})
            or (ex.seq_region_end > {1} and ex.seq_region_end < {2}))
            and tr.transcript_id = extr.transcript_id
            and extr.exon_id = ex.exon_id
            and tr.biotype != 'retained_intron';
        '''.format(chr, region_start, region_end)

        cursor = db.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute(exon_region_query)
        return cursor

def get_ccds_genes():
    ccds_genes_query = \
    '''
    select distinct
         gsid.stable_id as 'gene',
         seq_region.name as 'chr',
         gene.seq_region_start as 'start',
         gene.seq_region_end as 'end'
         seq_region.seq_region_id,
         seq_region.seq_region_strand as 'strand'
    from
         transcript as ts,
         object_xref as objxr,
         xref as xr,
         gene_stable_id as gsid,
         gene,
         seq_region
    where
         xr.external_db_id = 3800 and
         objxr.xref_id = xr.xref_id and
         objxr.ensembl_id = ts.transcript_id and
         objxr.ensembl_object_type = "Transcript" and
         gsid.gene_id = ts.gene_id and
         gene.gene_id = ts.gene_id and
         seq_region.seq_region_id = gene.seq_region_id
    order by gene.seq_region_id, gene.seq_region_start, gene.seq_region_end;
    '''

    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    cursor.execute(ccds_genes_query)
    return cursor

def get_ccds_exons(max_size):
    ccds_exons_query = \
    '''
    select distinct
         seq_region.name as 'chr',
         ex.seq_region_start as 'start',
         ex.seq_region_end as 'end',
         ex.stable_id as 'exon',
         GROUP_CONCAT(xr.display_label) as 'CCDS',
         GROUP_CONCAT(exts.rank) as 'ranks',
         ex.seq_region_end - ex.seq_region_start + 1 as 'len',
         ex.seq_region_strand as 'strand'
    from
         transcript as ts,
         exon as ex,
         exon_transcript as exts,
         xref as xr,
         object_xref as objxr,
         seq_region
    where
         ts.transcript_id = objxr.ensembl_id and     #transcript to extref
         objxr.ensembl_object_type = "Transcript" and
         objxr.xref_id = xr.xref_id and              #extref to extref type
         xr.external_db_id = 3800 and                #ref type is CCDS
         exts.transcript_id = ts.transcript_id and   #transcript to exon
         ex.phase = 0 and                            #phase
         ex.end_phase = 0 and
         ((ex.seq_region_end - ex.seq_region_start + 1) > 12 AND
          (ex.seq_region_end - ex.seq_region_start + 1) < %d) and
         ex.exon_id = exts.exon_id  and              #exon details
         seq_region.seq_region_id = ex.seq_region_id #seq region name
    GROUP BY ex.stable_id
    LIMIT 100;
    ''' % (max_size)

    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    cursor.execute(ccds_exons_query)
    return cursor

def main():
    print "Running through all exons...\n"
    categorize_exons()

if __name__ == "__main__":
    main()
