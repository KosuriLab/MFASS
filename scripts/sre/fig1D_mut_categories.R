load_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
  if(length(new_pkgs)) install.packages(new_pkgs)
  for(pkg in pkgs){
    suppressWarnings(suppressMessages(library(pkg, character.only = T)))
  }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'gridExtra', 'grid')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.png'

###############################################################################
# Read in data
###############################################################################
data <- read.table('../../processed_data/sre/sre_data_clean.txt', 
                   sep = '\t', header = T, 
                   colClasses = c('sub_id' = 'character')) %>%
  filter(rep_quality == 'high')

# read in re-scored reference file
updated_ref <- read.csv('../../ref/sre/sre_ref_rescored.txt',
                        header = T, sep = '\t',
                        colClasses = c('sub_id' = 'character'))

data <- data %>% 
  left_join(dplyr::select(updated_ref, id, exon_seq:correct_don_score), by = 'id')

###############################################################################
# Mutation categories
###############################################################################

esr_categories <- c('rmv_Ke2011_ESE', 'clst_Ke2011_ESE', 'rmv_Ke2011_ESS', 
                    'clst_Ke2011_ESS')
esr_labels <- c('weaken ESEs', 'destroy strongest ESE', 'weaken ESSs', 
                'destroy strongest ESS')

splice_site_categories <- c( 'weak_spl_a', 'weak_spl_d', 'p_weak_spl', 
                             'no_spl_a', 'no_spl_d',
                             'same_splice_a', 'same_splice_d',
                             'rmv_me_splice_acceptor', 'rmv_me_splice_donor', 
                             'csplice_a', 'csplice_d')
splice_site_labels <- c('weaken acceptor', 'weaken donor', 
                        'weaken donor + acceptor', 
                        'destroy acceptor', 'destroy donor',
                        'same score acceptor', 'same score donor',
                        'weaken spurious acceptor', 'weaken spurious donor',
                        'destroy spurious acceptor', 'destroy spurious donor')

intron_categories <- c('clst_Vlkr07_AICS', 'clst_Vlkr07_DICS', 
                       'rmv_Vlkr07_AICS', 'rmv_Vlkr07_DICS')
intron_labels <- c('weaken intronic conserved (acceptor side)', 
                   'weaken intronic conserved (donor side)', 
                   'destroy intronic conserved (acceptor side)', 
                   'destroy intronic conserved (donor side)')

random_exon_categories <- c('rnd_exon_1nt', 'rnd_exon_2nt', 
                            'rnd_exon_3nt', 'rnd_exon_5nt', 
                            'aggr_exon')
random_exon_labels <- c('random 1nt exon', 'random 2nt exon', 
                        'random 3nt exon', 'random 5nt exon',
                        'aggressive exon (only syn. mut.)')

random_intron_categories <- c('rnd_intron_1nt', 'rnd_intron_2nt', 
                              'rnd_intron_3nt', 'rnd_intron_5nt', 
                              'aggr_intron', 'p_aggr_intr', 'aggr_both')
random_intron_labels <- c('random 1nt intron', 'random 2nt intron', 
                          'random 3nt intron', 'random 5nt intron',
                          'aggressive intron', 'aggr. + random intron', 
                          'aggr. intron + exon')

other_categories <- c('cnsrv_1nt', 'cnsrv_3nt', 
                      'RBPmats', 'variation')
other_labels <- c('conserved 1nt', 'conservered 3nt', 
                  'destroy RBP motifs', 'dbSNPs')

gg <- data %>% 
  filter(nat_index_smn1 >= 0.50, !is.na(category)) %>%
  mutate(sdv_smn1 = ifelse(dpsi_smn1 <= -0.50, T, F)) %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>%
  group_by(category_fctr) %>% 
  summarise(num_sdv_smn1 = length(which(sdv_smn1 == T)),
            category_num_smn1 = n()) %>% 
  mutate(percent_sdv_smn1 = num_sdv_smn1/category_num_smn1 * 100) %>%
  ggplot(aes(x = category_fctr, y = percent_sdv_smn1)) +
  geom_col(alpha = 1, width = 0.8, fill = "blue") +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 105) +
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 8, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = "% SDV")

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/', 
              'fig1d_mut_categories', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = 100)

