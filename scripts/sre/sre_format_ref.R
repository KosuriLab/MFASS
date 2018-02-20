# Format SRE reference into distinct fields and convert genomic coordinates
# to lastest version

load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

###############################################################################
# File formatting
# read in ref, format ID field
ref <- read.table('../../ref/sre/sre_ref.txt', sep = '\t',
                  col.names = c('header', 'seq')) %>%
    # rowwise() %>%
    # mutate(id = as.character(unlist(strsplit(header, split = ' '))[[1]])) %>%
    # distinct(id, .keep_all = T) %>%
    # small substitutions to make separating easier
    mutate(header = gsub('cat= ', 'cat=', header))

# format nat and mut ref separately
nat_ref <- ref %>%
    filter(!grepl('mut', header)) %>%
    separate(header, into = c('id', 'chr', 'strand', 'length', 'ccds'), sep = ' ') %>%
    mutate(seq_type = 'nat')

mut_ref <- ref %>%
    filter(grepl('mut', header)) %>%
    # don't modify equal sign for control sequences because they don't have any 
    # information for those fields
    mutate(header = ifelse(grepl('CTRL', header), header, gsub('= ', '=', header)),
           header = gsub(', ', ',', header),
           header = gsub(';', '', header)) %>%
    separate(header, c('id', 'seq_type', 'chr', 'strand', 'length', 'ccds', 'category',
                       'loc', 'str', 'mutations', 'num_changes', 'orig_score',
                       'new_score', 'criteria_met', 'max_iter_reached',
                       'no_new_mutants'), sep = ' ') %>% 
    filter(!grepl('_000', id))

# more reference formatting
ref <- bind_rows(nat_ref, mut_ref) %>%
    distinct(id, .keep_all = T) %>% 
    separate(chr, c('chr', 'region'), sep = ':') %>%
    separate(region, c('start', 'end'), sep = '-', fill = 'right', convert = T) %>%
    separate(id, c('ensembl_id', 'sub_id'), sep = '_', remove = F) %>%
    # if natural sequence and no sub id, add 000
    mutate(sub_id = ifelse(is.na(sub_id), '000', sub_id)) %>%
    # get rid of leftover field identifiers
    mutate(strand = gsub('strand=', '', strand),
           length = gsub('len=', '', length),
           chr = gsub('chr', '', chr),
           ccds = gsub('ccds=', '', ccds),
           category = gsub('cat=', '', category),
           loc = gsub('loc=', '', loc),
           mutations = gsub('set=', '', mutations),
           num_changes = as.numeric(gsub('changes=', '', num_changes)),
           orig_score = as.numeric(gsub('orig_score=', '', orig_score)),
           new_score = as.numeric(gsub('new_score=', '', new_score)),
           criteria_met = as.logical(gsub('criteria_met=', '', criteria_met)),
           max_iter_reached = as.logical(gsub('MaxIterReached=', '', max_iter_reached)),
           no_new_mutants = as.logical(gsub('NoNewMutants=', '', no_new_mutants))) %>%
    # separate mutant range into separate columns
    separate(loc, c('mut_start', 'mut_end'), sep = ':', remove = F, convert = T) %>%
    # update seq type for controls
    mutate(seq_type = ifelse(grepl('CTRL', ensembl_id), 'control', seq_type),
           seq_type = ifelse(grepl('RAND', ensembl_id), 'control', seq_type)) %>%
    # finally, separate length into intron-exon-intron lengths
    extract(length, c("intron1_len","exon_len","intron2_len"),
            "([[:alnum:]]+).([[:alnum:]]+).([[:alnum:]]+)", convert = T) %>%
    arrange(ensembl_id, sub_id)

write.table(ref, '../../ref/sre/sre_ref_formatted.txt', sep = '\t',
row.names = F, col.names = F)

###############################################################################
# Convert genomic coordinates
###############################################################################

system('bash ./sre_ref_liftover.sh')

###############################################################################
# Join converted coordinates
###############################################################################
print('Converting genomic coordinates...')
ref <- ref %>%
    left_join(read.table('../../processed_data/sre/sre_ref_liftover.bed', 
                         sep = '\t', header = F, 
                         col.names = c('chr', 'start_hg38', 'end_hg38', 'id'))  %>%
                  select(-chr), by = 'id') %>% 
    mutate(id = ifelse(grepl('_', id), id, paste(ensembl_id, sub_id, sep = '_')))

                 
write.table(ref, '../../ref/sre/sre_ref_formatted_converted.txt',
  sep = '\t', row.names = F, col.names = T)

###############################################################################
# Update Ensembl exon IDs
###############################################################################
print('Updating Ensembl exon IDs...')
# ref <- read.table('../../ref/sre/sre_ref_formatted_converted.txt',
#                   sep = '\t', header = T)

# release 89 5/17 based on hg38
ensembl <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
attributes <- c('ensembl_gene_id', 'description', 'chromosome_name', 
                'start_position', 'end_position', 'strand',
                'ensembl_transcript_id', 'transcript_start', 'transcript_end', 
                'ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end',
                'is_constitutive', 'rank', 'phase', 'end_phase')

# grab all exons
all_exon_ids <- biomaRt::getBM(attributes = attributes, mart = ensembl)

# use start and end coordinates to find updated ID
exon_coords <- ref %>% 
    filter(seq_type == 'nat') %>% 
    distinct(ensembl_id, .keep_all = T) %>% 
    mutate(exon_start_hg38 = start_hg38 + intron1_len,
           exon_end_hg38 = exon_start_hg38 + exon_len - 1) %>% 
    select(exon_id_old = ensembl_id, exon_start_hg38, exon_end_hg38)

exon_coords_update <- exon_coords %>% 
    left_join(all_exon_ids, by = c('exon_start_hg38' = 'exon_chrom_start', 
                                   'exon_end_hg38' = 'exon_chrom_end')) %>% 
    rename(new_exon_id = ensembl_exon_id, gene_start_hg38 = start_position, 
           gene_end_hg38 = end_position)

# we only tested in-phase exons, keep these
exon_coords_update_inphase <- exon_coords_update %>% 
    filter(phase == 0, end_phase == 0)

write.table(exon_coords_update_inphase, '../../ref/exon_ids_updated.txt', 
            sep='\t', row.names = F)
