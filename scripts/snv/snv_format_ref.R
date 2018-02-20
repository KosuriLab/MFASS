# Format smv reference into distinct fields and convert genomic coordinates
# to lastest version

load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'Biostrings', 'stringr')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

###############################################################################
# File formatting
###############################################################################
ref <- read.table('../../ref/snv/snv_ref_all.txt', sep = '\t', 
                  col.names = c('header', 'sequence')) %>% 
    # remove spaces after equal signs and colons
    mutate(header = gsub(': ', ':', gsub('= ', '=', header))) %>% 
    # separate fields
    separate(header, c('id', 'region', 'strand', 'length', 'ref_allele', 'alt_allele', 
                       'snp_position', 'vcf_id', 'rel_position'), sep = ' ') %>% 
    separate(region, c('chr', 'region'), sep = ':') %>%
    separate(region, c('start', 'end'), sep = '-', fill = 'right', convert = T) %>%
    separate(id, c('ensembl_id', 'sub_id'), sep = '_', remove = F) %>%
    # get rid of leftover field identifiers
    mutate(strand = gsub('strand=', '', strand),
           length = gsub('len=', '', length),
           ref_allele = gsub('ref:', '', gsub('ref=', '', ref_allele)),
           alt_allele = gsub('alt=', '', alt_allele),
           # make sure this is numeric
           start = as.numeric(start),
           snp_position = as.numeric(gsub('pos=', '', snp_position)),
           vcf_id = as.numeric(gsub('vcf-id=', '', vcf_id)), 
           rel_position = as.numeric(gsub('rel_pos=', '', rel_position))) %>%
    # finally, separate length into intron-exon-intron lengths
    extract(length, c("intron1_len","exon_len","intron2_len"),
            "([[:alnum:]]+).([[:alnum:]]+).([[:alnum:]]+)", convert = T)

ctrl <- ref %>% 
    # the BRK controls have different formatting for reference alleles, 
    # separate these and deal with separately
    filter(endsWith(id, 'BRK')) %>% 
    # controls do not have SNP positions or VCF ids but have an additional set of reference and 
    # alternate alleles which occupy these columns. Combine this information into the reference 
    # and alternate allele columns
    mutate(snp_position = gsub(';', '', gsub('ref=', '', snp_position)),
           ref_allele = paste0(ref_allele, snp_position),
           snp_position = NA,
           vcf_id = gsub('alt=', '', vcf_id),
           alt_allele = paste0(alt_allele, ';', vcf_id),
           vcf_id = NA) %>% 
    # combine back with reference
    bind_rows(filter(ref, !endsWith(id, 'BRK'))) %>% 
    # now re-format and convert the rest of the necessary columns
    mutate(snp_position = as.numeric(gsub('pos=', '', snp_position)),
           rel_position = as.numeric(gsub('rel_pos=', '', rel_position)),
           vcf_id = as.numeric(gsub('vcf-id=', '', vcf_id)))

write.table(ref, '../../ref/snv/snv_ref_formatted.txt', sep = '\t',
            row.names = F, col.names = F, quote = F)

###############################################################################
# convert genome coordinates
###############################################################################
system('bash ./snv_ref_liftover.sh')

# join
ref <- ref %>% 
    left_join(read.table('../../processed_data/snv/snv_ref_liftover.bed', sep = '\t',
                         col.names = c('chr', 'start_hg38_0based', 'end_hg38_0based', 'id')) 
              %>% select(-chr), by = 'id') %>% 
    left_join(read.table('../../processed_data/snv/snv_snp_liftover.bed', sep = '\t',
                         col.names = c('chr', 'snp_position_hg38_0based_start', 
                                       'snp_position_hg38_1based', 'id')) %>% 
                  select(-chr) %>% 
                  mutate(snp_position_hg38_0based_end = snp_position_hg38_1based), 
              by = 'id')


# get relative position of SNP mutations and where they fall (intron/exon)
rel_position <- function(intron1_len, exon_len, intron2_len, strand, rel_position) {
    in_interval <- function(start, end, x) {
        if (start <= x & x <= end) { return(TRUE) }
        else {return(FALSE)}
    }
    if (strand == '-' | strand == -1) {
        # set appropriate lengths
        upstr_intron_len <- intron2_len
        downstr_intron_len <- intron1_len
    }
    else{
        upstr_intron_len <- intron1_len
        downstr_intron_len <- intron2_len
    }
    
    regions <- data.frame(label = c('upstr_intron', 'exon', 'downstr_intron'),
                          start = c(1, upstr_intron_len + 1, upstr_intron_len + exon_len + 1),
                          end = c(upstr_intron_len, upstr_intron_len + exon_len, 
                                  upstr_intron_len + exon_len + downstr_intron_len))
    # select region the SNP falls in
    label <- regions %>% 
        rowwise() %>% 
        filter(in_interval(start, end, x = rel_position))
    
    start <- label$start[1]
    end <- label$end[1]
    label <- label$label[1]
    
    # get distance from end of upstream intron/exon boundary
    boundary <- upstr_intron_len
    distance <- rel_position - boundary
    
    # normalize to feature length
    if (label == 'downstr_intron') {
        scaled_distance <- 1 + (distance - exon_len) / (end - start)
    }
    else {
        scaled_distance <- distance / (end - start + 1)
    }
    
    # additionally, find distance from respective intron/exon boundary
    # left side of boundary is negative, right side is positive
    if (label == 'upstr_intron') {
        rel_pos_feature <- rel_position - end - 1
    }
    if (label == 'exon') {
        # closer to downstream intron
        if ( rel_position - start + 1 >= end - rel_position + 1) {
            rel_pos_feature <- rel_position - end - 1 # negative position, left side of boundary
        }
        else {
            # closer to upstream intron, right side of boundary, positive
            rel_pos_feature <- rel_position - start + 1
        }
    }
    if (label == 'downstr_intron') {
        rel_pos_feature <- rel_position - start + 1
    }
    return(paste(label, rel_pos_feature, scaled_distance, sep = ':'))
}

ref <- ref %>% 
    mutate(rel_position = ifelse(strand == '+', snp_position - start + 1, 
                                 end - snp_position + 1))

ref <- ref %>% 
    select(intron1_len, exon_len, intron2_len, strand, rel_position, id) %>%
    na.omit() %>%
    rowwise() %>%
    mutate(rel_position_info = rel_position(intron1_len, exon_len, 
                                            intron2_len, strand, rel_position)) %>%
    separate(rel_position_info, 
             c('label', 'rel_position_feature', 'rel_position_scaled'), 
             sep = ':', convert = T) %>%
    select(id, label, rel_position_feature, rel_position_scaled) %>%
    left_join(ref, ., by = 'id') %>% distinct()

write.table(ref, '../../ref/snv/snv_ref_formatted_converted.txt', sep = '\t',
            quote = F, row.names = F)

###############################################################################
# It is harder to synthesize A's, so oligos in the reference file are the strand
# with the lower A count. Let's read in the original reference file without 
# flipped sequence
###############################################################################

ref <- read.table('../../ref/snv/snv_ref_formatted_converted.txt', 
                  sep = '\t', header = T)
snv_ref_original_seq <- read.table('../../ref/snv/snv_ref_original_seq.txt',
                                    sep = '\t', header = T) %>% 
  select(id, original_seq = sequence)

ref <- ref %>%
  left_join(snv_ref_original_seq, by = 'id') %>%
  mutate(original_seq = ifelse(is.na(original_seq), sequence, original_seq),
         mixed_seq = sequence) %>%
  arrange(id, sub_id)

# make column for natural sequence
ref <- ref %>% 
    group_by(ensembl_id) %>%
    mutate(natural_seq = ifelse(any(sub_id == '000'), 
                                original_seq[sub_id == '000'],
                                NA)) %>%
    ungroup()

write.table(ref, '../../ref/snv/snv_ref_formatted_converted_original_seq.txt', 
            sep = '\t', quote = F, row.names = F)
