### Add external functional data ###

load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.png'

# read in data
data <- read.table('../../processed_data/snv/snv_data_clean.txt',
                   sep = '\t', header = T)

# remove irrelevant columns
data <- data %>% 
    select(-mixed_seq, -nat_seq)

###############################################################################
# extract conservation for each SNP
###############################################################################
# output bed file of SNP positions to get conservation
write.table(file = '../../processed_data/snv/snp_positions_hg38.bed', 
            x = data %>% 
                select(chr, snp_position_hg38_0based_start, 
                       snp_position_hg38_0based_end, id) %>% 
                na.omit(), 
            sep = '\t', col.names = F, row.names = F, quote = F)

# phastCons
system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/snv/snp_positions_hg38.bed',
             '../../processed_data/snv/snp_position_cons_scores.bed'))

# read in phastCons
phastCons <- read.table('../../processed_data/snv/snp_position_cons_scores.bed', 
                        sep = '\t', header = F,
                        col.names = c('name', 'size', 'bases_covered', 
                                      'snp_sum', 'mean0', 'mean_phastCons_score')) %>% 
    filter(bases_covered != 0)

data <- data %>% 
  left_join(select(phastCons, id = name, mean_phastCons_score), by = 'id') 


# run phyloP
system(paste('bash',
             '../run_phyloP.sh',
             '../../processed_data/snv/snp_positions_hg19.bed',
             '../../processed_data/snv/snp_position_phyloP_scores.bed'))
# read in phyloP
phylop <- read.table('../../processed_data/snv/snp_position_phyloP_scores.bed', 
                        sep = '\t', header = F,
                        col.names = c('name', 'size', 'bases_covered', 
                                      'snp_sum', 'mean0', 'phylop_score')) %>% 
    filter(bases_covered != 0)

data <- data %>% 
    left_join(select(phylop, id = name, phylop_score), by = 'id') 


###############################################################################
# grab ExAC annotation for all sequences
###############################################################################
ref <- read.table('../../ref/snv/snv_ref_formatted_converted_original_seq.txt',
                  sep = '\t', header = T)
# write BED file of SNP positions
write.table(ref %>%
                filter(!is.na(snp_position)) %>% 
                mutate(chr = gsub('chr', '', chr),
                       # ExAC built on hg19, supply 1-based region so tabix works properly
                       snp_region = paste0(chr, ':', snp_position, '-',  snp_position + 1)) %>% 
                select(snp_region) ,
            file = '../../processed_data/snv/tabix_input_snp_regions.txt',
            quote = F, row.names = F, col.names = F, sep = '\t')

# write SNP file for 0-based hg19 coordinates
write.table(file = '../../processed_data/snv/snp_positions_hg19.bed', 
            x = data %>% 
                mutate(snp_position_0based_start = snp_position - 1) %>% 
                select(chr, snp_position_0based_start, snp_position, id) %>% 
                na.omit(),
            sep = '\t', col.names = F, row.names = F, quote = F)

# this takes awhile
# system(paste('while read line; do tabix',
#              '../../ref/snv/ExAC.r0.3.1.sites.vep.vcf.gz',
#              '$line >> ../../ref/snv/snp_snv_annot_all.txt;',
#              'done < ../../processed_data/snv/tabix_input_snp_regions.txt'))

# read in annotation and parse
exac_annot <- read.table('../../ref/snv/snp_snv_annot_all.txt', 
                         sep = '\t', header = F)
colnames(exac_annot) <- c('chr', 'position', 'snp_id', 'ref_allele', 
                          'alt_allele', 'quality', 'filter', 'info')

# grab allele frequency column separately
exac_annot <- exac_annot %>% 
    separate(info, into = c('other', 'AF'), sep = 'AF=', remove = F) %>% 
    select(-other) %>% 
    separate(AF, into = c('AF', 'other'), sep = ';', convert = T) %>% 
    select(-other)

exac_annot <- exac_annot %>% 
    # only care about CSQ field, split by it and get rid of everything before
    separate(info, into = c('other', 'CSQ'), sep = 'CSQ=') %>% 
    select(-other) %>% 
    # get rid of all fields after CSQ
    separate(CSQ, into = c('CSQ', 'other'), sep = ';') %>% 
    select(-other)

# parse CSQ column, consequence annotations from Ensembl VEP (Variant Effect 
# Predictor). Description of fields taken from header of vcf file
csq_fields <- 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF|context|ancestral'
csq_fields <- unlist(strsplit(csq_fields, split = '[|]'))
exac_annot <- separate(exac_annot, CSQ, into = csq_fields, sep = '[|]')

exac_annot <- exac_annot %>% 
    # give each alternate allele its own entry
    mutate(alt_allele = strsplit(alt_allele, split = ','),
           AF = strsplit(AF, split = ',')) %>%
    distinct(chr, position, ref_allele, alt_allele, Consequence, .keep_all = T) %>%
    # make sure alt allele matches allele the annotation is describing
    filter(alt_allele == Allele) %>% 
    mutate(AF = as.numeric(AF)) %>%
    # rearrange
    select(chr:Allele, annot = Consequence, IMPACT:AF)

# join exac data
data <- data %>%
    mutate(chr = gsub('chr', '', chr)) %>% 
    left_join(select(exac_annot, alt_allele = Allele, chr, position, snp_id, annot,
                     impact = IMPACT, symbol = SYMBOL, symbol_source = SYMBOL_SOURCE, HGNC_ID,
                     ensembl_gene_id = Gene, ensembl_feature_id = Feature,
                     exon_number = EXON, intron_number = INTRON, Existing_variation,
                     allele_num = ALLELE_NUM, SIFT, PolyPhen, GMAF:EA_MAF, LoF,
                     context, ancestral, AF),
              by = c('chr', 'snp_position' = 'position', 'alt_allele')) %>% 
    distinct(id, annot, .keep_all = T)

###############################################################################
# re-score variants using Ensembl VEP
###############################################################################
# write VCF file for Ensembl VEP (1-based)
data %>% 
    filter(category == 'mutant') %>% 
    select(chr, snp_position_hg38_1based, id, ref_allele, alt_allele) %>% 
    write.table(file = '../../ref/snv/snv.vcf', 
                sep='\t', quote = F, row.names = F, col.names = F)

# this takes awhile
# system(paste('../../bin/ensembl-vep/vep -i ../../ref/snv/snv.vcf',
#              '--cache --dir_cache ../../ref/vep/'))
# system('mv variant_effect_output.txt ../../processed_data/snv/snv_vep.txt')
# system('rm -f variant_effect_output.txt_summary.html')

vep <- read.table('../../processed_data/snv/snv_vep.txt', sep = '\t', comment.char = '#')
colnames(vep) <- c('id', 'location', 'alt_allele', 'ensembl_gene_id_vep', 
                   'ensembl_transcript_id', 'feature_type', 'consequence', 
                   'cDNA_position', 'CDS_position', 'protein_position', 
                   'amino_acids', 'codons', 'existing_variation', 'extra')

ensembl <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
gene_info <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'start_position', 'end_position',
                                           'external_gene_name', 'external_gene_source', 
                                           'gene_biotype', 'hgnc_id', 'hgnc_symbol'), 
                            mart = ensembl,
                            filters = 'ensembl_gene_id',
                            values = vep$ensembl_gene_id_vep)

gene_info <- gene_info %>% 
    mutate(gene_length = end_position - start_position)

vep <- vep %>% 
    left_join(select(gene_info, -start_position, -end_position),
              by = c('ensembl_gene_id_vep' = 'ensembl_gene_id'))

# filter out any entries with distance in extra column, we don't want distal transcripts
vep <- vep %>% 
    filter(!grepl('DISTANCE', extra)) %>% 
    separate(extra, into = c('impact', 'strand_vep'), sep = ';') %>% 
    select(-impact) %>% 
    mutate(strand_vep = gsub('STRAND=', '', strand_vep),
           strand_vep = gsub('-1', '-', strand_vep),
           strand_vep = gsub('1', '+', strand_vep))


data_with_duplicates <- left_join(select(data, id:HGNC_ID,
                                         ensembl_gene_id_exac = ensembl_gene_id,
                                         ensembl_feature_id:AF), 
                 select(vep, id, location, ensembl_gene_id_vep, external_gene_name,
                        external_gene_source, consequence, strand_vep,
                        gene_biotype, gene_length) %>% 
                     distinct(ensembl_gene_id_vep, location, .keep_all = T), 
                 by = 'id')

data_no_duplicates <- data_with_duplicates %>% 
    group_by(id) %>% 
    filter(n() == 1) %>% 
    ungroup()

duplicates <- data_with_duplicates %>% 
    group_by(id) %>% 
    filter(n() > 1)

# check if exac and vep ids agree, pick the matching one
data_no_duplicates <- duplicates %>%
    group_by(id) %>% 
    filter(ensembl_gene_id_exac == ensembl_gene_id_vep) %>% 
    bind_rows(data_no_duplicates, .)

# add back if VEP id present but no exac id
data_no_duplicates <- duplicates %>% 
    anti_join(., data_no_duplicates, by = 'id') %>% 
    filter(is.na(ensembl_gene_id_exac), !is.na(ensembl_gene_id_vep)) %>% 
    group_by(id) %>% 
    slice(1) %>% 
    ungroup() %>% 
    bind_rows(data_no_duplicates, .)

# check which ones are still missing
missing <- anti_join(data, data_no_duplicates, by = 'id') %>% 
    semi_join(data_with_duplicates, ., by = 'id') %>% 
    filter(strand == strand_vep)

# these missing ones have directly conflicting IDs, choose VEP, newer annotation
data_no_duplicates <- missing %>% 
    group_by(id) %>% 
    slice(1) %>% 
    ungroup() %>% 
    bind_rows(., data_no_duplicates)

# finally, check if VEP strand matches and fill in gene for wild-type sequence
data <- data_no_duplicates %>% 
    mutate(ensembl_gene_id = ifelse(strand == strand_vep, 
                                    ensembl_gene_id_vep, 
                                    ensembl_gene_id_exac)) %>% 
    select(-ensembl_gene_id_exac, -ensembl_gene_id_vep, -strand_vep) %>% 
    group_by(ensembl_id) %>% 
    mutate(symbol = ifelse(is.na(symbol), symbol[category == 'mutant'], symbol),
           ensembl_gene_id = ifelse(is.na(ensembl_gene_id), 
                                    ensembl_gene_id[category == 'mutant'],
                                    ensembl_gene_id),
           external_gene_name = ifelse(is.na(external_gene_name), 
                                       unique(external_gene_name), 
                                       external_gene_name)) %>% 
    ungroup()

# fill in external_gene_name for those that are missing but have other entries 
# with ensembl_gene_id with information
data <- data %>% 
    group_by(ensembl_gene_id) %>% 
    filter(any(!is.na(external_gene_name))) %>% 
    mutate(external_gene_name = ifelse(is.na(external_gene_name), 
                                       unique(external_gene_name), 
                                       external_gene_name)) %>% 
    ungroup()
    

rm(data_with_duplicates, duplicates, data_no_duplicates, missing)

# join updated Ensembl exon IDs
updated_ids <- read.table('../../ref/exon_ids_updated.txt',
                          sep = '\t', header = T)

# exons can fall in more than one transcript, annotations are dependent on both
# transcript and exon ID
exon_data <- data %>%
    select(id, chr, snp_position_hg38_1based, ref_allele, alt_allele,
           exon_id_old = ensembl_id) %>%
    left_join(select(updated_ids, exon_id_old, exon_id_new = new_exon_id,
                     ensembl_gene_id, ensembl_transcript_id, is_constitutive,
                     rank), by = c('exon_id_old'))

exon_data <- exon_data %>%
    left_join(vep, by = c('id', 'alt_allele'))
#     distinct()

# let's try to filter it so there is one transcript per exon. For now, let's 
# just choose the transcript where the consequence is the most severe
# convert to factor, order levels according to severity can be found here:
# (http://www.ensembl.org/info/genome/variation/predicted_data.html)
consequence_levels <- c('transcript_ablation', 'splice_acceptor_variant', 
                        'splice_donor_variant', 'stop_gained', 'frameshift_variant', 
                        'stop_lost', 'start_lost', 'transcript_amplification', 
                        'inframe_insertion', 'inframe_deletion', 'missense_variant', 
                        'protein_altering_variant', 'splice_region_variant',
                        'incomplete_terminal_codon_variant', 'stop_retained_variant',
                        'synonymous_variant', 'coding_sequence_variant', 
                        'mature_miRNA_variant', '5_prime_UTR_variant', 
                        '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 
                        'intron_variant', 'NMD_transcript_variant',
                        'non_coding_transcript_variant', 'upstream_gene_variant', 
                        'downstream_gene_variant', 'TFBS_ablation', 
                        'TFBS_amplification', 'TF_binding_site_variant',
                        'regulatory_region_ablation', 'feature_elongation', 
                        'regulatory_region_variant', 'feature_truncation', 
                        'intergenic_variant')

# some entries have multiple consequences, create new row for each
exon_data <- exon_data %>% 
    mutate(consequence = strsplit(consequence, ',')) %>% 
    unnest(consequence) %>% 
    mutate(consequence = factor(consequence,
                                levels = consequence_levels)) %>% 
    group_by(id) %>% 
    arrange(consequence) %>%
    filter(row_number() == 1) %>% 
    distinct(id, consequence, .keep_all = T) %>% 
    ungroup() %>% 
    arrange(id)

# join to original data
data <- data %>%
    select(-consequence) %>% 
    left_join(select(exon_data, id, alt_allele, exon_id_new, consequence), 
              by = c('id', 'alt_allele')) 

# sdv_by_category <- data %>% 
#     group_by(consequence) %>% 
#     summarise(num_sdv = length(which(strong_lof == T)),
#               category_num = n()) %>% 
#     mutate(percent_sdv = num_sdv / category_num) %>% 
#     arrange(desc(num_sdv))

###############################################################################
# CADD annotation
###############################################################################
# submitted VCF file (snv.vcf) at http://cadd.gs.washington.edu/score to 
# score variants
cadd_annot <- read.table('../../processed_data/snv/snp_cadd_annot.txt', 
                         sep = '\t', header = F, 
                         col.names = c('chr', 'snp_position', 'ref_allele', 
                                       'alt_allele', 'cadd_score', 'phred')) %>%
    distinct()

data <- data %>%
    left_join(select(cadd_annot, -ref_allele, -phred), 
              by = c('chr', 'snp_position', 'alt_allele'))

###############################################################################
# FATHMM-MKL annotation
###############################################################################
# http://fathmm.biocompute.org.uk/fathmmMKL.html
# based on hg19
data %>% 
    filter(category == 'mutant') %>%
    select(chr, snp_position, ref_allele, alt_allele) %>% 
    write.table('../../processed_data/snv/fathmm_input.txt',
                sep = '\t', quote = F)

# input file submitted on web server
fathmm <- read.table('../../processed_data/snv/snv_fathmm_scores.tsv',
                     sep = '\t', header = F, skip = 1,
                     col.names = c('chr', 'snp_position', 'ref_allele',
                                   'alt_allele', 'fathmm_noncoding_score',
                                   'noncoding_group', 'fathmm_coding_score',
                                   'coding_group', 'warning'))

data <- data %>% 
    left_join(select(fathmm, -noncoding_group, -coding_group, -warning), 
              by = c('chr', 'snp_position', 'ref_allele', 'alt_allele'))

###############################################################################
# fitCons annotation
###############################################################################
# hg19 fitCons, 0-based bed file
#http://compgen.cshl.edu/fitCons/0downloads/tracks/i6/scores/
system(paste('bigWigAverageOverBed ../../ref/fc-i6-0.bw',
             '../../processed_data/snv/snp_positions_hg19.bed',
             '../../processed_data/snv/snp_fitCons_scores.bed'))

fitcons <- read.table('../../processed_data/snv/snp_fitCons_scores.bed',
                      sep = '\t', header = F,
                      col.names = c('name', 'size', 'bases_covered', 
                                    'snp_sum', 'mean0', 'fitCons_score'))
data <- data %>% 
    left_join(select(fitcons, id = name, fitCons_score), by = 'id')

###############################################################################
# DANN annotation
###############################################################################
#https://cbcl.ics.uci.edu/public_data/DANN/. # whole file is 102G
# system(paste('while read line; do tabix ../../ref/DANN_whole_genome_SNVs.tsv.bgz $line',
#              '>> ../../processed_data/snv/snp_dann_annot.txt done',
#              '< ../../processed_data/snv/tabix_input_snp_regions.txt'))

dann <- read.table('../../processed_data/snv/snp_dann_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
    left_join(dann, by = c('chr', 'snp_position', 'ref_allele', 'alt_allele'))

###############################################################################
# LINSIGHT annotation
###############################################################################
# http://compgen.cshl.edu/~yihuang/LINSIGHT/
# hg19/hg19
system(paste('bigWigAverageOverBed ../../ref/linsight.bw',
             '../../processed_data/snv/snp_positions_hg19.bed',
             '../../processed_data/snv/snp_linsight_score.bed'))

linsight <- read.table('../../processed_data/snv/snp_linsight_score.bed',
                       sep = '\t', header = F,
                       col.names = c('id', 'size', 'bases_covered', 
                                     'snp_sum', 'mean0', 'linsight_score'))

data <- data %>% 
    left_join(select(linsight, id, linsight_score), by = 'id')

data <- data %>% 
    distinct(.keep_all = T)

###############################################################################
# probability of gene being loss-of-function intolerant
###############################################################################

pli_table <- read.table('../../ref/snv/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt',
                        sep = '\t', header = T)

data <-  data %>%
    left_join(select(pli_table, gene, pLI), by = c('external_gene_name' = 'gene')) %>%
    mutate(pLI_high = ifelse(pLI >= 0.90, TRUE, FALSE))

# tmp <- data %>% 
#     left_join(select(pli_table, gene, pLI), by = c('symbol' = 'gene')) %>% 
#     mutate(pLI_high = ifelse(pLI >= 0.90, TRUE, FALSE))
# 
# tmp2 <-  data %>% 
#     left_join(select(pli_table, gene, pLI), by = c('external_gene_name' = 'gene')) %>% 
#     mutate(pLI_high = ifelse(pLI >= 0.90, TRUE, FALSE))
# 
# tmp3 <- left_join(select(tmp, id, pLI_high_symbol = pLI_high, symbol),
#                   select(tmp2, id, pLI_high_name = pLI_high, external_gene_name))
# 
# table(select(tmp3, pLI_high_symbol, pLI_high_name))

write.table(data, '../../processed_data/snv/snv_func_annot.txt',
            sep = '\t', quote = F, row.names = F)

# write simpler data frame of variant list
data %>% 
    filter(category != 'control') %>% 
    dplyr::select(internal_id = id, category, ref_allele, alt_allele, 
                  start, end, strand,
                  intron1_len, exon_len, intron2_len,
                  current_ensembl_exon_id = exon_id_new, 
                  ensembl_gene_id, ensembl_transcript_id,
                  chr, strand, start_hg19_0based = start, end_hg19_0based = end,  
                  start_hg38_0based, end_hg38_0based, 
                  snp_position_hg19_0based = snp_position,
                  snp_position_hg38_0based_start, snp_position_hg38_0based_end,
                  label, seq = original_seq, 
                  index = v2_index, nat_index = nat_v2_index, 
                  dpsi = v2_dpsi, is_sdv = strong_lof) %>% 
    write.table('../../processed_data/snv/snv_list.txt', 
                sep = '\t', row.names = F, quote = F)


tested_intron <- c('ENSE00001264396_002', 'ENSE00001320758_023', 
                   'ENSE00001111517_006', 'ENSE00001125865_010')
tested_syn <- c('ENSE00000904134_014', 'ENSE00001320758_011',
                'ENSE00001367617_012', 'ENSE00002386775_008')
tested_sa <- c('ENSE00001635327_019', 'ENSE00001117772_021',
               'ENSE00001649206_010', 'ENSE00001007663_004')
tested_sd <- c('ENSE00001033448_030', 'ENSE00001137258_018',
               'ENSE00001144033_005', 'ENSE00001649206_024', 'ENSE00001639529_009')
tested_splice_region <- c('ENSE00000815598_010', 'ENSE00001093636_009',
                          'ENSE00001152991_007', 'ENSE00001093061_012' ,
                          'ENSE00001125865_034') 
tested_missense <- c('ENSE00000815598_013', 'ENSE00001732647_010', 
                     'ENSE00001735670_006', 'ENSE00001308061_007') 

tested_exons <- data.frame(id = c(tested_intron, tested_syn,
                                  tested_sa, tested_sd,
                                  tested_splice_region, tested_missense),
                           type = c(rep('intron' , length(tested_intron)),
                                    rep('synonymous', length(tested_syn)),
                                    rep('splice_acceptor', length(tested_sa)),
                                    rep('splice_donor', length(tested_sd)),
                                    rep('splice_region', length(tested_splice_region)),
                                    rep('missense', length(tested_missense))))
tested_exons <- tested_exons %>% 
    left_join(select(data, id, ensembl_gene_id, symbol, gene_length, 
                     snp_position, snp_position_hg38_0based_start, exon), by = 'id') %>% 
    arrange(gene_length)

write.table(sdv %>% 
                mutate(snp_position_0based_start = snp_position - 1, 
                       chrom = paste0('chr', chr)) %>% 
                select(chrom, snp_position_0based_start, snp_position, id),
            '../../processed_data/snv/sdv_hg19.bed', sep = '\t', quote = F, 
            row.names = F, col.names = F)

write.table(sdv %>% 
                mutate(chrom = paste0('chr', chr)) %>% 
                select(chrom, snp_position_hg38_0based_start, snp_position_hg38_0based_end, id),
            '../../processed_data/snv/sdv_hg38.bed', sep = '\t', quote = F, 
            row.names = F, col.names = F)



