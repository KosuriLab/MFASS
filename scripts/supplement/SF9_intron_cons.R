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

options(stringsAsFactors = F, warn = -1, warnings = -1)
plot_format <- '.png'
plot_format_main <- '.tiff'

data <- read.table('../../processed_data/snv/snv_data_clean.txt', 
                   sep = '\t', header = T)

###############################################################################
# calculate genomic coordinates for the upstream and downstream introns for
# conservation
###############################################################################
intron_coords <- data %>%
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index>=0.5) %>%
    distinct(ensembl_id, .keep_all = T) %>%
    mutate(intron1_start = start_hg38_0based,
           intron1_end = start_hg38_0based + intron1_len - 1,
           intron2_start = start_hg38_0based + intron1_len + exon_len,
           intron2_end = end_hg38_0based) %>%
    dplyr::select(id, ensembl_id, chr, strand, intron1_len, intron2_len, intron1_start:intron2_end)

intron_coords %>%
    mutate(upstr_intron = ifelse(strand == '+',
                                 paste(intron1_start, intron1_end, sep = '-'),
                                 paste(intron2_start, intron2_end, sep = '-')),
           downstr_intron = ifelse(strand == '+',
                                   paste(intron2_start, intron2_end, sep = '-'),
                                   paste(intron1_start, intron1_end, sep = '-'))) %>%
    dplyr::select(chr, upstr_intron, downstr_intron, ensembl_id) %>%
    gather('feature_type', 'region', upstr_intron:downstr_intron) %>%
    separate(region, into = c('start', 'end'), sep = '-') %>%
    mutate(id = paste(ensembl_id, feature_type, sep = '-')) %>%
    dplyr::select(chr, start, end, id) %>%
    write.table(file = '../../processed_data/snv/snv_intron_coords.bed',
                sep = '\t', quote = F, row.names = F, col.names = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/snv/snv_intron_coords.bed',
             '../../processed_data/snv/snv_intron_cons_scores_all.bed'))

intron_cons <- read.table('../../processed_data/snv/snv_intron_cons_scores_all.bed', 
                          sep = '\t', header = F,
                          col.names = c('name', 'size', 'bases_covered', 
                                        'snp_sum', 'mean0', 'mean_cons_score')) %>% 
    filter(bases_covered != 0) %>% 
    separate(name, into = c('ensembl_id', 'feature_type'), sep = '-') %>% 
    mutate(feature_type = paste0(feature_type, '_mean_cons')) %>% 
    spread(feature_type, mean_cons_score)

data <- data %>% 
    left_join(dplyr::select(intron_cons, ensembl_id, upstr_intron_mean_cons) %>% 
                  na.omit(),
              by = 'ensembl_id') %>% 
    left_join(dplyr::select(intron_cons, ensembl_id, downstr_intron_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id')

###############################################################################
# Let's get coordinates for 100bp of the flanking intron instead of shorter 
# version we synthesized.
###############################################################################
intron_coords %>%
    mutate(intron1_start = intron1_start - (100 - intron1_len + 1),
           intron2_end = intron2_end + (100 - intron2_len + 1),
           upstr_intron_100 = ifelse(strand == '+',
                                     paste(intron1_start, intron1_end, sep = '-'),
                                     paste(intron2_start, intron2_end, sep = '-')),
           downstr_intron_100 = ifelse(strand == '+',
                                       paste(intron2_start, intron2_end, sep = '-'),
                                       paste(intron1_start, intron1_end, sep = '-'))) %>%
    dplyr::select(chr, upstr_intron_100, downstr_intron_100, ensembl_id) %>%
    gather('feature_type', 'region', upstr_intron_100:downstr_intron_100) %>%
    separate(region, into = c('start', 'end'), sep = '-') %>%
    mutate(id = paste(ensembl_id, feature_type, sep = '-')) %>%
    dplyr::select(chr, start, end, id) %>%
    write.table(file = '../../processed_data/snv/snv_intron_coords_100.bed',
                sep = '\t', quote = F, row.names = F, col.names = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/snv/snv_intron_coords_100.bed',
             '../../processed_data/snv/snv_intron_cons_scores_100_all.bed'))

intron_cons_100 <- read.table('../../processed_data/snv/snv_intron_cons_scores_100_all.bed', 
                              sep = '\t', header = F,
                              col.names = c('name', 'size', 'bases_covered', 
                                            'snp_sum', 'mean0', 'mean_cons_score')) %>% 
    filter(bases_covered != 0) %>% 
    separate(name, into = c('ensembl_id', 'feature_type'), sep = '-') %>% 
    mutate(feature_type = paste0(feature_type, '_mean_cons')) %>% 
    spread(feature_type, mean_cons_score)

data <- data %>% 
    left_join(dplyr::select(intron_cons_100, ensembl_id, upstr_intron_100_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id') %>% 
    left_join(dplyr::select(intron_cons_100, ensembl_id, downstr_intron_100_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id')

data <- data %>% 
    mutate(upstr_intron_len = ifelse(strand == '+', intron1_len, intron2_len),
           downstr_intron_len = ifelse(strand == '+', intron2_len, intron1_len)) %>%
    filter(!is.na(v2_dpsi), category == "mutant", nat_v2_index >= 0.5) %>%
    mutate(`Excluded intronic context (up to 100bp, acceptor)` = 
             upstr_intron_mean_cons - upstr_intron_100_mean_cons,
           `Excluded intronic context (up to 100bp, donor)` = 
             downstr_intron_mean_cons - downstr_intron_100_mean_cons) %>%
    mutate(`Included intronic context (acceptor)` = 
             upstr_intron_mean_cons,
           `Included intronic context (donor)` = 
             downstr_intron_mean_cons)

data2 <- data %>%
    gather(key = 'acceptor_conservation_type', value = 'acceptor_phastCons', 
           c(`Excluded intronic context (up to 100bp, acceptor)`, 
             `Included intronic context (acceptor)`)) %>%
    gather(key = 'donor_conservation_type', value = 'donor_phastCons', 
           c(`Excluded intronic context (up to 100bp, donor)`, 
             `Included intronic context (donor)`)) 

# supplement, synthetic short intron conservation vs. 100 bp conservation
# upstream
gg <- ggplot(data2, aes(acceptor_phastCons)) +
  geom_density(alpha = 0.5, 
               aes(fill = acceptor_conservation_type)) +
  labs(x = 'average phastCons score\nsplice acceptor side'
  ) +
  scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(
    legend.position = 'none',
    axis.title.x = element_text(size = 20, vjust = -2),
    axis.title.y = element_text(size = 20, vjust = +4),
    axis.text.x = element_text(size = 14, color = 'grey20'),
    axis.text.y = element_text(size = 14, color = 'grey20'),
    axis.ticks.x = element_line(color = 'grey50'),
    axis.ticks.y = element_line(color = 'grey50'),
    axis.line.x = element_line(color = 'grey50'),
    axis.line.y = element_line(color = 'grey50'),
    plot.margin = unit(c(2,2,3,3),"mm"))
gg

ggsave(paste0('../../figs/supplement/SF9A_acceptor_conservation', plot_format),
       gg, width = 5, height = 5, dpi = 300, scale = 1.3)

# downstream
gg <- ggplot(data2, aes(donor_phastCons)) +
  geom_density(alpha = 0.5, 
               aes(fill = donor_conservation_type)) +
  labs(x = 'average phastCons score\nsplice donor side') +
  scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(
    legend.position = 'none',
    axis.title.x = element_text(size = 20, vjust = -2),
    axis.title.y = element_text(size = 20, vjust = +4),
    axis.text.x = element_text(size = 14, color = 'grey20'),
    axis.text.y = element_text(size = 14, color = 'grey20'),
    axis.ticks.x = element_line(color = 'grey50'),
    axis.ticks.y = element_line(color = 'grey50'),
    axis.line.x = element_line(color = 'grey50'),
    axis.line.y = element_line(color = 'grey50'),
    plot.margin = unit(c(2,2,3,3),"mm"))
gg

ggsave(paste0('../../figs/supplement/SF9B_snv_donor_conservation', plot_format),
       gg, width = 5, height = 5, dpi = 300, scale = 1.3)

########################
# Legend for Figure
########################

# grab legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

gg <- ggplot(data2, aes(acceptor_phastCons)) +
  geom_density(alpha = 0.5, 
               aes(
                 # color = acceptor_conservation_type,
                 fill = acceptor_conservation_type)) +
  # geom_abline(intercept = 0, slope = 1, type = 'dashed') +
  labs(x = 'average phastCons score\nsplice acceptor side'
  ) +
  scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  # scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20, vjust = -2),
    axis.title.y = element_text(size = 20, vjust = +4),
    axis.text.x = element_text(size = 14, color = 'grey20'),
    axis.text.y = element_text(size = 14, color = 'grey20'),
    axis.ticks.x = element_line(color = 'grey50'),
    axis.ticks.y = element_line(color = 'grey50'),
    axis.line.x = element_line(color = 'grey50'),
    axis.line.y = element_line(color = 'grey50'),
    plot.margin = unit(c(2,2,3,3),"mm"))
gg

legend <- g_legend(gg)
tiff(paste0('../../figs/supplement/SF9_legend_acceptor', plot_format_main), 
     width = 100, height = 14, units = 'mm', res = 300)
grid.newpage()
grid.draw(legend)
dev.off()

gg <- ggplot(data2, aes(donor_phastCons)) +
  geom_density(alpha = 0.5, 
               aes(fill = donor_conservation_type)) +
  labs(x = 'average phastCons score\nsplice donor side') +
  scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 20, vjust = -2),
    axis.title.y = element_text(size = 20, vjust = +4),
    axis.text.x = element_text(size = 14, color = 'grey20'),
    axis.text.y = element_text(size = 14, color = 'grey20'),
    axis.ticks.x = element_line(color = 'grey50'),
    axis.ticks.y = element_line(color = 'grey50'),
    axis.line.x = element_line(color = 'grey50'),
    axis.line.y = element_line(color = 'grey50'),
    plot.margin = unit(c(2,2,3,3),"mm"))
gg

legend <- g_legend(gg)
tiff(paste0('../../figs/supplement/SF9_legend_donor', plot_format_main), 
     width = 100, height = 14, units = 'mm', res = 300)
grid.newpage()
grid.draw(legend)
dev.off()
