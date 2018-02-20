import sys
import copy

import Bio.Alphabet
from Bio import SeqIO

import cfg
import feature
import score
import mutate
#from feature import SMSeqRecord

def main(seq_file,options):

    #sequence file (fasta)
    sys.stderr.write("Reading file: "+seq_file+"\n")

    try:
        seq_record = SeqIO.parse(seq_file, "gb").next()
    except ValueError:
        #SeqBuilder can screw up spacing...
        feature.fix_sdb_genbank(seq_file)
        seq_record = SeqIO.parse(seq_file, "gb").next()

    sys.stderr.write("  Parsing record: "+seq_record.id+"\n")

    #convert to an SMSeqRecord, which has extra functions that we need...
    #import pdb; pdb.set_trace()
    #seq_record.__class__ = SMSeqRecord
    #seq_record.populate_attribs()

    #APE can screw up alphabets...
    seq_record.seq.alphabet = Bio.Alphabet.generic_dna

    seq_record = feature.remove_annotations(seq_record)
    seq_record = score.annotate_splice_signals(seq_record)
    seq_score = score.score_seq_record(seq_record)

    min_score = seq_score
    best_record = seq_record


    if options['enable_mutate']:

        #start from scratch 'mutate_record_max' times, and within those
        #seperate tries, try to remove motif sets iteratively by conservative
        #single base mutations in individual motifs that reduce motif scores.
        for seq_iter in range(cfg.mutate_record_max):

            sys.stderr.write("Iteration: {0}\tBest: {1}\tOrig:{2}\n".format(
                seq_iter,min_score,seq_score))

            if min_score == 0:
                break

            edited_sr = copy.deepcopy(seq_record)
            last_kept_sr = copy.deepcopy(seq_record)

            #iterate motif removal 'mutate_record_iter' times. At the end
            #of each iteration, don't use the motif removal step if it has
            #increased the total record score (sum of all motif scores)
            for i in range(cfg.mutate_record_iter):
                last_kept_sr = copy.deepcopy(edited_sr)

                if min_score == 0:
                    break

                edited_sr = mutate.remove_all_motifs(edited_sr)
                edited_sr = feature.remove_annotations(edited_sr)
                edited_sr = score.annotate_splice_signals(edited_sr)

                #import pdb; pdb.set_trace()

                if (score.score_seq_record(edited_sr) <
                    score.score_seq_record(last_kept_sr)):
                    last_kept_sr = edited_sr
                    sys.stderr.write("  New edit score: {0} was kept.\n"
                        .format(score.score_seq_record(edited_sr)))
                else:
                    sys.stderr.write(
                        "  New edit score: {0} was higher than {1}.\n"
                        .format(
                                score.score_seq_record(edited_sr),
                                score.score_seq_record(last_kept_sr)))

            #check final total record score
            #edit_score = scoreSeqRecord(edited_sr)

            if score.score_seq_record(last_kept_sr) < min_score:
                min_score = score.score_seq_record(last_kept_sr)
                best_record = last_kept_sr
            else:
                del edited_sr

    #as a final check, do a whole-sequence translation to ensure no amino acids
    #have changed...
    mutate.check_all_translations(best_record,seq_record)

    if options['output_fasta']:
        print best_record.format('fasta')
    else:
        print best_record.format('gb')


#    ###############################
#    # Get splice element motifs
#    ###############################
#
#    #find and print short splicing motifs
#    motifSets = load_motif_sets()
#    foundMotifs = find_strings(seq_record,motifSets)
#    for mSet in foundMotifs:
#        for string in foundMotifs[mSet]:
#            for pos in foundMotifs[mSet][string]:
#                print "\t".join([mSet+":"+string,
#                                 str(pos[0]),
#                                 str(pos[1]),
#                                 '-'])
#
#    ###############################
#    # etcs
#    ###############################
#
#    #index potential donor splice sites
#    donorLocations = where(MEdonorScores > 0)
#    for donorLoc in sorted(donorLocations[0]):
#        #import pdb; pdb.set_trace()
#        print "\t".join(["donor    ",
#            str(donorLoc+maxEntDonor[0]),
#            str(donorLoc+maxEntDonor[1]),
#            str(MEdonorScores[donorLoc])])
#
#    #index potential acceptor splice sites
#    acceptorLocations = where(MEacceptorScores > 0)
#    for acceptorLoc in sorted(acceptorLocations[0]):
#
#        print "\t".join(["acceptor",
#            str(acceptorLoc+maxEntAcceptor[0]),
#            str(acceptorLoc+maxEntAcceptor[1]),
#            str(MEacceptorScores[acceptorLoc])])
#
#        #find branch sites that are applicable:
#        bsSearchStart = acceptorLoc - cfg.bsUpstreamSearch
#        if bsSearchStart < 0: bsSearchStart = 0
#        bsSearchStr = str(seq_record.seq[bsSearchStart:acceptorLoc])
#        bsList = call_find_branch_site(bsSearchStr)
#
#        for bs in bsList:
#            print "\t".join(["  branch_site",
#                str(bs['start']+bsSearchStart),
#                str(bs['end']+bsSearchStart),
#                "mut["+str(bs['score'])+"]"])
#
#            pptSearchRegion = ((bs['end']+bsSearchStart),acceptorLoc)
#            ppt = find_ppt(seq_record.seq[slice(*pptSearchRegion)])
#            if len(ppt) > 0:
#                print "\t".join(["  ppt     ",
#                    str(ppt['start']+pptSearchRegion[0]),
#                    str(ppt['end']+pptSearchRegion[0]),
#                    "len["+str(ppt['end']-ppt['start']+1)+"]"])
sys.stderr.write("Done.\n")

if __name__ == "__main__":
        #####program options
    seq_file = sys.argv[1]
    options = {}

    #enable mutation of cryptic sites so that they will be removed
    options['enable_mutate'] = sys.argv.count('-mutate') > 0
    options['output_fasta'] = sys.argv.count('-fasta') > 0

    main(seq_file,options)





