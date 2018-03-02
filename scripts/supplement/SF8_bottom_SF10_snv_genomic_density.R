## This script generates supplementary figure 8 (bottom panel) and 10 
## for SNV genome wide density 

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

plot_format <- '.tiff'


###############################################################################
# SNV density, genome-wide
###############################################################################
# extract intronic regions from Gencode annotation file
system('bash ./extract_intronic_regions.sh')

# grabs all SNVs from ExAC, convert from vcf to bed, uses intersectBed against
# bed file for all introns/exons to find which genomic feature each SNV falls in
# this takes a few minutes
# system('bash ./find_snv_genomic_feature.sh')

nat_snps <- read.table('../../ref/snv/snv_genomic_overlap.bed.gz',
                       sep = '\t', header = F)
colnames(nat_snps) <-  c('snp_chrom', 'snp_start', 'snp_end', 'snp_ID', 
                         'snp_score', 'strand', 'chrom', 'start', 'end', 'id')

# let's get rid of the score column and any SNVs without a strand
nat_snps <- nat_snps %>% 
    filter(strand != '.') %>% 
    select(-snp_score)

# calculate absolute position and scaled position. For introns, let's only 
# consider sequences n bp (e.g. 100) from start, and assign
# the first nbp as downstream intron and last nbp as upstream intron, to simplify 
# graphing
nat_snps <- nat_snps %>% 
    mutate(rel_position = ifelse(strand == '1', 
                                 snp_start - start + 1, 
                                 end - snp_start + 1),
           feature_len = end - start + 1)

n <- 100
nat_snps <- nat_snps %>% 
    mutate(label = case_when(.$id != 'intron' ~ 'exon',
                             .$id == 'intron' & .$rel_position <= n ~ 'downstr_intron',
                             .$id == 'intron' & .$feature_len - .$rel_position <= n ~ 'upstr_intron',
                             TRUE ~ 'intron'))


### relative position of SNVs, genome-wide ###
nat_snps <- nat_snps %>% 
    mutate(rel_position = ifelse(label == 'upstr_intron', 
                                 # make upstream intron negative for graphing purposes
                                 feature_len - rel_position, 
                                 rel_position),
           label_len = case_when(.$label == 'upstr_intron' | .$label == 'downstr_intron' ~ n,
                                 TRUE ~ .$feature_len))

nat_snps <- nat_snps %>% 
    mutate(rel_position_scaled = rel_position / label_len,
           rel_position_scaled = ifelse(label == 'upstr_intron', 
                                        -1 * rel_position_scaled, 
                                        rel_position_scaled),
           rel_position_scaled = ifelse(label == 'downstr_intron', 
                                        1 + rel_position_scaled, 
                                        rel_position_scaled))

nat_snp_summary <- nat_snps %>% 
    mutate(rel_pos_binned = cut(rel_position_scaled, breaks = seq(-1, 2, 0.01), include.lowest = T)) %>% 
    filter(label != 'intron') %>% 
    group_by(rel_pos_binned) %>% 
    summarise(num_snps = n())

nat_snp_summary %>% 
    ggplot(aes(rel_pos_binned, 0.5)) + geom_tile(aes(fill = log(num_snps))) + 
    viridis::scale_fill_viridis() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          legend.text = element_text(size=10),
          legend.title = element_text(size=8)) +
    labs(x = 'relative position', y = '', fill = 'log average\nnumber of\nSNPs',
         title = 'All exonic SNVs, intronic SNVs within 100bp') 

ggsave(paste0('../../figs/supplement/SF8bottom_genome_SNVs_rel_position', plot_format),
       height = 1.5, width = 8, unit = 'in')


### average number of SNVs intron vs. exon, genome-wide ###
nat_snps %>% 
    # combine upstream and downstream into one for graph purposes
    mutate(label = ifelse(grepl('intron', label), 'intron', 'exon'),
           rel_position_scaled = rel_position / feature_len,
           rel_pos_binned = cut(rel_position_scaled, 
                                breaks = seq(0, 1, 0.01), 
                                include.lowest = T)) %>% 
    group_by(rel_pos_binned, label) %>% 
    summarise(num_snps = n()) %>% 
    ggplot(aes(rel_pos_binned, num_snps)) + geom_point() +
    facet_grid(. ~ label) +
    coord_cartesian(ylim = c(0,38500)) +
    scale_y_continuous(breaks = seq(10000, 30000, by = 10000), expand = c(0,0)) +
    labs(x = 'scaled position', y = 'average\nnumber of SNVs') +
    theme(
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        legend.position = 'none',
        axis.title.x = element_text(size = 14, vjust = -2), 
        axis.title.y = element_text(size = 14, vjust = +4),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = 'grey20'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'grey50'),
        axis.line.x = element_line(color = 'grey50'),
        axis.line.y = element_line(color = 'grey50'),
        plot.margin = unit(c(2,2,3,3),"mm")) 

ggsave(paste0('../../figs/supplement/SF10_genome_num_SNVs_exon_vs_intron', plot_format),
       height = 3, width = 4, unit = 'in')
