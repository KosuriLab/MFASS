import Bio
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import re
import sys
import copy

import UserString

#import splicemod utils
sys.path.append("/Users/dbgoodman/Dropbox/intron/code/splicemod_dg/src")
import cfg
import feature
import score
#from feature import SMSeqRecord
import mutate
import util


#fasta_path = \
#  "/Users/dbgoodman/Dropbox/intron/data/candidate_exons/ccds_CNE_final.fasta"

#out_path = \
#  "/Users/dbgoodman/Dropbox/intron/data/candidate_exons"

# def score_exons(fasta_file, output_path):

# def mutate_exons(chr_ranges,):

def main(options):

    seq_count = -1

    if options['print_all_scores']:
        output_file = "%s/%s.txt" % (options['output_path'],
                                     'all_candidate_scores.txt')
        output_handle = open(output_file, "w")

    records = [record for record in \
               Bio.SeqIO.parse(options['input_file'],"fasta")]


    for record in records:

        seq_count += 1

        cut_sites = {'gcgatcgc':'AsiSI','cctgcagg':'SbfI','ggcgcgcc':'AscI',
                     'ttaattaa':'PacI'}

        #add the description to the record annotations:
        record.annotations.update(map(lambda x: x.split('='),
                                      record.description.split(' ')[1:]))

        if not options['print_all_scores'] and \
            record.annotations["range"] not in options['chr_ranges']:
            continue

        for site_seq in cut_sites.keys():
            if site_seq.lower() in str(record.seq.lower()):
                print 'WARNING: cut site %s found in this record!' % cut_sites[site_seq]


        #parse the values in the name of each record
        #example: hg19_ct_UserTrack_3545_CCDS46737.1_8/23_60_78_40
        #separated by underscores, they are:
        #              name, exon rank, upstream, exon, downstream
        annot_list = ['gene_id','exon_rank','upstr_intron_size', \
                     'exon_size','downstr_intron_size']
        for k,v in zip(annot_list,record.name.split('_')[-5:]):
            record.annotations[k] = v

        #convert numbers to int()
        for k,v in record.annotations.iteritems():
            if v.isdigit(): record.annotations[k] = int(v)

        #set the gene id & name to the CCDS id + the exon rank
        record.id = record.annotations["gene_id"]
        record.name = record.id

        #set the record alphabet
        record.seq.alphabet = Bio.Alphabet.DNAAlphabet()

        #create the sequence features:

        #calculate the feature start/ends
        upstr_loc = (0,record.annotations["upstr_intron_size"])
        exon_loc = (upstr_loc[1],upstr_loc[1]+record.annotations["exon_size"])
        downstr_loc = (exon_loc[1],exon_loc[1]+
                       record.annotations["downstr_intron_size"])

        #upstream intron
        record.features.append( \
            Bio.SeqFeature.SeqFeature(
                Bio.SeqFeature.FeatureLocation(*upstr_loc),
                type= "intron",
                strand= record.annotations["strand"] == '+'
            )
        )

        #exon
        record.features.append( \
            Bio.SeqFeature.SeqFeature(
                Bio.SeqFeature.FeatureLocation(*exon_loc),
                type= "exon",
                strand= record.annotations["strand"] == '+'
            )
        )

        #downstream intron
        record.features.append( \
            Bio.SeqFeature.SeqFeature(
                Bio.SeqFeature.FeatureLocation(*downstr_loc),
                type= "intron",
                strand= record.annotations["strand"] == '+'
            )
        )

        #convert to an SMSeqRecord, which has extra functions that we need...
        #record.__class__ = SMSeqRecord
        record.populate_attribs()
        record = score.annotate_splice_signals(record)
        record.get_correct_motifs()

        #populate cryptic motif list
        record.cryptic_motifs = []
        for feat in record.features:
            if not "evidence" in feat.qualifiers: continue
            if "function" in feat.qualifiers and \
                feat.qualifiers["function"].count('correct splice signal') > 0:
                continue
            record.cryptic_motifs.append(feat)

        #print all cryptic and correct scores to stdout
        if options['print_all_scores']:

            exon_length_string = str.join('-',map( str, \
                    [record.annotations["upstr_intron_size"],
                    record.annotations["exon_size"],
                    record.annotations["downstr_intron_size"]]))

            correct_motif_score_sum = \
                sum(map(lambda x: float(x.qualifiers["evidence"]),
                        record.correct_motifs))

            cryptic_motif_score_sum = \
                sum(map(lambda x: float(x.qualifiers["evidence"]),
                        record.cryptic_motifs))

            rec_str = \
                  "size= %s\trange= %s\tstrand= %s\tcorrect motifs= %d\t"+ \
                  "correct score= %2.2f\tcryptic motifs= %d\t"+ \
                  "cryptic score=%2.2f" % ( \
                    exon_length_string,
                    record.annotations["range"],
                    record.annotations["strand"],
                    len(record.correct_motifs),
                    correct_motif_score_sum,
                    len(record.cryptic_motifs),
                    cryptic_motif_score_sum)

            output_handle.write(rec_str+"\n")

        #if no -printall flag, then instead create genbank files of ranges
        # and add motifs to them
        else:
            output_file = "%s/%s.fasta" % (options['output_path'],
                                           record.annotations["gene_id"])
            output_handle = open(output_file, "w")

            #TODO: probably want to make these seqs an input via a separate
            #genbank file, so that we can use different cut sites on a dime
#            record.seq = 'ATTGCGATCGC' + record.seq
#            record.seq = record.seq + 'CCTGCAGGAAT'
#
#            #add features that correspond to the cut sites
#            record.features.append(
#                Bio.SeqFeature.SeqFeature(
#                    Bio.SeqFeature.FeatureLocation(3,11),
#                    type= "misc_feature",
#                    qualifiers = {'note':'AsiSI Site'},
#                    strand= record.annotations["strand"] == '+'
#                )
#            )
#
#            record.features.append(
#                Bio.SeqFeature.SeqFeature(
#                    Bio.SeqFeature.FeatureLocation(190,198),
#                    qualifiers = {'note':'SbfI Site'},
#                    type= "misc_feature",
#                    strand= record.annotations["strand"] == '+'
#                )
#            )

            #finally, write the output file
            SeqIO.write(record, output_handle, "genbank")
            output_handle.close()

            #testing motif mutation:
            record.correct_motifs.sort( \
                key=lambda seq_feat: seq_feat.extract_pos())

            iranges = cfg.gen_mut_ranges

            #TODO: I'd like to consolidate all the maxent data into a class...

            #generate mutants in iranges and add them to the 'mutants' qualifier
            #in each of the correct motif features (acceptor and donor)
            for feat in record.correct_motifs:
                try:
                    invariant_mask = cfg.maxEntInvariants[
                                                    feat.qualifiers['label'][0]]
                except:
                    raise ValueError('''Masking failed; is the feature's label
                                        qualifier correct?''')

                mutate.generate_feature_mutants(record,feat,
                                                iranges,invariant_mask)


            #pick pairs of feature mutants (closest to middles of iranges) and
            #generate and print new records with those mutants


            don_pos = slice(*record.correct_motifs[0].extract_pos())
            don_oseq = str(record.seq[don_pos]).lower()
            acc_pos = slice(*record.correct_motifs[1].extract_pos())
            acc_oseq = str(record.seq[acc_pos]).lower()

            #include first and last codon that are in splice signals
            exon = record.seq[don_pos.stop-3:acc_pos.start+3]

            translation = Seq(str(exon)).translate().tostring()

            print "########\n#Finished Mutant: \n>seq%d-O name=%s acc=%s don=%s trans=%s" % \
                (seq_count,record.description,
                 float(record.correct_motifs[0].qualifiers['evidence']),
                 float(record.correct_motifs[1].qualifiers['evidence']),
                 translation)

            pseq = str(record.seq).lower()

            def seq_with_spaces(seq,slice1,slice2):
                return seq[:slice1[0].start].lower() + ' ' \
                       + slice1[1][:-3] + ' ' + slice1[1][-3:] + ' ' \
                       + seq[slice1[0].stop:slice2[0].start].lower() + ' ' \
                       + slice2[1][:3] + ' '  + slice2[1][3:] + ' ' \
                       + seq[slice2[0].stop:].lower()

            print "%s" % seq_with_spaces(str(record.seq),
                                         (don_pos,don_oseq),
                                         (acc_pos,acc_oseq))
            #import pdb; pdb.set_trace()

            print "# Originals:"
            print "# \t\t%s\t%s\t(%2.2f)\t(%2.2f)" % \
                (don_oseq, acc_oseq,
                 float(record.correct_motifs[0].qualifiers['evidence']),
                 float(record.correct_motifs[1].qualifiers['evidence'])
                )

            rng_names = ('W','L','H','S') #weak, low, high, strong

            for irange in iranges:
                print "\n# For score range %s:" % str(irange)
                don = record.correct_motifs[0]
                acc = record.correct_motifs[1]
                don_mut_tuples = don.qualifiers['mutants'][irange].items()
                acc_mut_tuples = acc.qualifiers['mutants'][irange].items()

                #calculate the distance to the center of the specified range
                # for each of the donor and acceptor tuples
                calc_range_dist = lambda dt: util.range_dist(dt[1],irange)

                #TODO: Put num_range_members in cfg
                rng_mbrs = 2

                #sort the donor and acceptors by distance from the center of the
                # score range
                don = sorted(don_mut_tuples, key=calc_range_dist)[:(rng_mbrs)]
                acc = sorted(acc_mut_tuples, key=calc_range_dist)[:(rng_mbrs)]

                #take one from each of these sorted sets to make pairs
                for pair in range(rng_mbrs):
                    print "#  Pair %d:\t%s\t%s\t(%2.2f)\t(%2.2f)" \
                    % (pair,don[pair][0],acc[pair][0],don[pair][1],acc[pair][1])

                    pseq = str(record.seq).lower()

                    print ">seq%d-%s-%d name=%s pair=%d irange=%s acc=%s don=%s trans=%s" % \
                        (seq_count, rng_names[iranges.index(irange)], pair,
                         record.description, pair, irange,
                         don[pair][1],acc[pair][1],
                         translation)

                    print "%s\n" % seq_with_spaces(str(record.seq),
                                         (don_pos,don[pair][0]),
                                         (acc_pos,acc[pair][0]))


    if options['print_all_scores']:
        output_handle.close()


# ============
# = END MAIN =
# ============

if __name__ == "__main__":
    #####program options
    options = {}

    options['print_all_scores'] = sys.argv.count('-printall') > 0
    options['input_file'] = sys.argv[sys.argv.index('-i')+1]
    options['output_path'] = sys.argv[sys.argv.index('-o')+1]
    if not options['print_all_scores']:
        options['chr_ranges'] = \
            sys.argv[sys.argv.index('-ranges')+1].split(',')

    main(options)


