### This script makes supplementary figures for genome-wide statistics ###

load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1, scipen = 10000)
plot_format <- '.tiff'

data <- read.table('../../processed_data/snv/snv_data_clean.txt', 
                   sep = '\t', header = T)

###############################################################################
# Conservation of in-frame vs. out-of-frame exons (whole genome)
###############################################################################
# release 89 5/17 based on hg38
ensembl <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
attributes <- c('ensembl_gene_id', 'description', 'chromosome_name', 
                'start_position', 'end_position', 'strand', 'ensembl_transcript_id', 
                'transcript_start', 'transcript_end', 'ensembl_exon_id', 
                'exon_chrom_start', 'exon_chrom_end', 'is_constitutive', 'rank', 
                'phase', 'end_phase')

all_exon_ids <- biomaRt::getBM(attributes = attributes, mart = ensembl)

# calculate various coordinates necessary to determine relative scaled position
# of SNVs within intron/exon
all_exon_ids <- all_exon_ids %>% 
    mutate(intron1_end = exon_chrom_start - 1,
           intron1_start = intron1_end - 99,
           intron2_start = exon_chrom_end + 1,
           intron2_end = intron2_start + 99,
           upstr_intron_start = ifelse(strand == 1, intron1_start, intron2_start),
           upstr_intron_end = ifelse(strand == 1, intron1_end, intron2_end),
           downstr_intron_start = ifelse(strand == 1, intron2_start, intron1_start),
           downstr_intron_end = ifelse(strand == 1, intron2_end, intron1_end),
           upstr_intron_len = upstr_intron_end - upstr_intron_start + 1,
           downstr_intron_len = downstr_intron_end - downstr_intron_start + 1,
           exon_len = exon_chrom_end - exon_chrom_start + 1)

inframe_outframe_exons <- all_exon_ids %>% 
    filter(phase != -1 & end_phase != -1) %>% 
    mutate(exon_type = case_when(.$phase == 0 & .$end_phase == 0 ~ 'inframe',
                                 TRUE ~ 'outframe')) %>% 
    distinct(ensembl_exon_id, .keep_all = T)

# let's get upstream intron coordinates first. We will generate a bed file
# containing position of each nucleotide in range so we can get conservation per
# nucleotide
inframe_outframe_exons %>% 
    filter(chromosome_name %in% c(1:22, 'X', 'Y')) %>%
    mutate(chromosome_name = paste0('chr', chromosome_name)) %>% 
    select(chromosome_name, upstr_intron_start, upstr_intron_end, ensembl_exon_id) %>% 
    group_by(ensembl_exon_id) %>% 
    mutate(position = list(seq(upstr_intron_start, upstr_intron_end, by = 1))) %>% 
    unnest(position) %>%
    ungroup() %>% 
    mutate(id = paste(ensembl_exon_id, position, sep = '_'),
           start = position,
           end = position + 1) %>% 
    select(chromosome_name, start, end, id) %>% 
    write.table(file = '../../processed_data/snv/nat_upstr_intron_positions.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/snv/nat_upstr_intron_positions.bed',
             '../../processed_data/snv/nat_upstr_intron_cons_scores_all.bed'))

# downstream intron
inframe_outframe_exons %>% 
    filter(chromosome_name %in% c(1:22, 'X', 'Y')) %>%
    mutate(chromosome_name = paste0('chr', chromosome_name)) %>% 
    select(chromosome_name, downstr_intron_start, downstr_intron_end, ensembl_exon_id) %>% 
    group_by(ensembl_exon_id) %>% 
    mutate(position = list(seq(downstr_intron_start,downstr_intron_end, by = 1))) %>% 
    unnest(position) %>%
    ungroup() %>% 
    mutate(id = paste(ensembl_exon_id, position, sep = '_'),
           start = position,
           end = position + 1) %>% 
    select(chromosome_name, start, end, id) %>% 
    write.table(file = '../../processed_data/snv/nat_downstr_intron_positions.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/snv/nat_downstr_intron_positions.bed',
             '../../processed_data/snv/nat_downstr_intron_cons_scores_all.bed'))

# exon
inframe_outframe_exons %>%                             
    filter(chromosome_name %in% c(1:22, 'X', 'Y')) %>%
    mutate(chromosome_name = paste0('chr', chromosome_name)) %>% 
    select(chromosome_name, exon_chrom_start, exon_chrom_end, ensembl_exon_id) %>% 
    group_by(ensembl_exon_id) %>% 
    mutate(position = list(seq(exon_chrom_start,exon_chrom_end, by = 1))) %>% 
    unnest(position) %>%
    ungroup() %>% 
    mutate(id = paste(ensembl_exon_id, position, sep = '_'),
           start = position,
           end = position + 1) %>% 
    select(chromosome_name, start, end, id) %>% 
    write.table(file = '../../processed_data/snv/nat_exon_positions.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/snv/nat_exon_positions.bed',
             '../../processed_data/snv/nat_exon_cons_scores_all.bed'))

# it takes awhile to read in the conservation score files and calculate summary,
# so we'll only do this if they don't already exist
if(!file.exists('../../processed_data/snv/nat_upstr_cons_summary.rds')){
    nat_upstr_cons <- read.table('../../processed_data/snv/nat_upstr_intron_cons_scores_all.bed', 
                                 sep = '\t', header = F,
                                 col.names = c('name', 'size', 'bases_covered', 
                                               'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>% 
        tidyr::separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = T) %>% 
        select(-(size:mean0)) %>% 
        left_join(select(inframe_outframe_exons, ensembl_exon_id, upstr_intron_len, 
                         upstr_intron_start, upstr_intron_end, strand, exon_type),
                  by = c('ensembl_id' = 'ensembl_exon_id')) %>%
        arrange(ensembl_id, position) %>% 
        mutate(position = as.numeric(position),
               rel_position = ifelse(strand == 1, 
                                     upstr_intron_end - position,
                                     position - upstr_intron_start),
               # upstream intron, keep them all negative,
               rel_position = -1 * rel_position,
               rel_position_scaled = rel_position / upstr_intron_len,
               rel_pos_binned = cut(rel_position_scaled, breaks = seq(-1, 0, 0.01))) %>% 
        group_by(rel_pos_binned, exon_type) %>%
        summarise(mean_cons_per_rel_pos = mean(mean_cons_score, na.rm = T))
    # save as RDS so levels in factor variable are saved correctly
    saveRDS(nat_upstr_cons, '../../processed_data/snv/nat_upstr_cons_summary.rds')
} else{
    nat_upstr_cons <- readRDS('../../processed_data/snv/nat_upstr_cons_summary.rds')
}

if(!file.exists('../../processed_data/snv/nat_downstr_cons_summary.rds')) {
    nat_downstr_cons <- read.table('../../processed_data/snv/nat_downstr_intron_cons_scores_all.bed', 
                                   sep = '\t', header = F,
                                   col.names = c('name', 'size', 'bases_covered', 
                                                 'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>%
        tidyr::separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = T) %>% 
        select(-(size:mean0)) %>% 
        left_join(select(inframe_outframe_exons, ensembl_exon_id, downstr_intron_len, 
                         downstr_intron_start, downstr_intron_end, strand, exon_type),
                  by = c('ensembl_id' = 'ensembl_exon_id')) %>%
        arrange(ensembl_id, position) %>% 
        mutate(position = as.numeric(position),
               rel_position = ifelse(strand == 1,
                                     position - downstr_intron_start,
                                     downstr_intron_end - position),
               rel_position_scaled = 1 + (rel_position / downstr_intron_len),
               rel_pos_binned = cut(rel_position_scaled, breaks = seq(1, 2, 0.01), 
                                    include.lowest = T)) %>% 
        group_by(rel_pos_binned, exon_type) %>%
        summarise(mean_cons_per_rel_pos = mean(mean_cons_score, na.rm = T))
    saveRDS(nat_downstr_cons, '../../processed_data/snv/nat_downstr_cons_summary.rds')
} else {
    nat_downstr_cons <- readRDS('../../processed_data/snv/nat_downstr_cons_summary.rds')
}

if(!file.exists('../../processed_data/snv/nat_exon_cons_summary.rds')) {
    nat_exon_cons <- read.table('../../processed_data/snv/nat_exon_cons_scores_all.bed', 
                                sep = '\t', header = F,
                                col.names = c('name', 'size', 'bases_covered', 
                                              'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>%
        tidyr::separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = T) %>% 
        select(-(size:mean0)) %>% 
        left_join(select(inframe_outframe_exons, ensembl_exon_id, exon_len, 
                         exon_chrom_start, exon_chrom_end, strand, exon_type),
                  by = c('ensembl_id' = 'ensembl_exon_id')) %>%
        arrange(ensembl_id, position) %>% 
        mutate(position = as.numeric(position),
               rel_position = ifelse(strand == 1,
                                     position - exon_chrom_start,
                                     exon_chrom_end - position),
               rel_position_scaled = rel_position / exon_len,
               rel_pos_binned = cut(rel_position_scaled, breaks = seq(0, 1, 0.01), 
                                    include.lowest = T)) %>% 
        group_by(rel_pos_binned, exon_type) %>%
        summarise(mean_cons_per_rel_pos = mean(mean_cons_score, na.rm = T))
    saveRDS(nat_exon_cons, '../../processed_data/snv/nat_exon_cons_summary.rds')
} else {
    nat_exon_cons <- readRDS('../../processed_data/snv/nat_exon_cons_summary.rds')
}

nat_cons <- bind_rows(nat_upstr_cons, nat_exon_cons, nat_downstr_cons)
nat_cons$rel_pos_binned <- factor(nat_cons$rel_pos_binned, levels = unique(nat_cons$rel_pos_binned))
nat_cons$exon_type <- factor(nat_cons$exon_type)
levels(nat_cons$exon_type) <- c('phase (0-0)', 'other phases')

# natural conservation between in-frame and out-of-frame exons, genome-wide
ggplot(nat_cons, aes(rel_pos_binned, mean_cons_per_rel_pos, color = exon_type)) +
    geom_point() + scale_color_manual(values = c('black', 'red')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(legend.position = c(0.85, 0.80)) +
    labs(x = '', 
         y = 'average\nphastCons score', color = 'exon phase') +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    coord_cartesian(ylim=c(0, 1)) +
    theme(legend.title = element_blank(),
        legend.position = c(0.675, 0.7),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = 'grey20'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'grey50'),
        axis.line.x = element_line(color = 'grey50'),
        axis.line.y = element_line(color = 'grey50'),
        legend.text = element_text(size = 12))

ggsave(paste0('../../figs/supplement/SF8_genome_exon_cons_inframe_vs_outframe', plot_format),
       height = 4, width = 5, unit = 'in')

save.image("../../processed_data/snv/snv_intron_cons.RData")
