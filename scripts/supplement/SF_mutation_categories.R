###############################################################################
# set-up
###############################################################################
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

plot_format_main <- '.tiff'
plot_format <- '.tiff'
hi_res <- 600

jitter_alpha <- 0.10

# custom color palette
source("../color_palette.R")
# specify new color palette
steps <- c("blue2", "cyan", "white", "yellow", "red2")
pal <- color.palette(steps, c(160, 1, 1, 160), space = "rgb")

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

# # calculate change in donor/acceptor score between mutant and wild-type 
# data <- data %>% 
#   group_by(ensembl_id) %>% 
#   filter(any(sub_id == '000')) %>% 
#   mutate(correct_acc_score_nat = correct_acc_score[sub_id == '000'],
#          correct_don_score_nat = correct_don_score[sub_id == '000'],
#          delta_acc_score = correct_acc_score - correct_acc_score_nat,
#          delta_don_score = correct_don_score - correct_don_score_nat,
#          # calculate fold change in score relative to wild-type
#          don_score_fold_change = (correct_don_score - correct_don_score_nat) / 
#            abs(correct_don_score_nat),
#          acc_score_fold_change = (correct_acc_score - correct_acc_score_nat) / 
#            abs(correct_acc_score_nat)) %>%
#   ungroup()

###############################################################################
# sre categories
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

# SMN1 intron backbone INCLUDED
gg <- data %>% 
  filter(seq_type == 'mut', nat_index_smn1 >= 0.5) %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
              # aes(color = nat_index_smn1)
              ) + 
  geom_boxplot(alpha = 0) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
                         # breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 10, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/smn1/', 
              'sre_smn1_all_categories_no_x_axis', plot_format_main), 
       gg_no_x_axis, width = 12, height = 3, dpi = hi_res)
ggsave(paste0('../../figs/sre/smn1/', 
              'sre_smn1_all_categories_no_x_axis_INCLUDED', plot_format_main), 
       gg_no_x_axis, width = 12, height = 3, dpi = hi_res)
ggsave(paste0('../../figs/sre/smn1/',
              'sre_smn1_all_categories', plot_format), 
       gg, width = 12, height = 5, dpi = hi_res)

# SMN1 intron backbone EXCLUDED
gg <- data %>% 
  filter(seq_type == 'mut', nat_index_smn1 < 0.5) %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
              # aes(color = nat_index_smn1)
  ) + 
  geom_boxplot(alpha = 0) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  # breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 10, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/smn1/', 
              'sre_smn1_all_categories_no_x_axis_EXCLUDED', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = 100)

# percent splice-disrupting mutants by category
sdv_by_category_smn1 <- data %>% 
  filter(nat_index_smn1 >= 0.75) %>%
  mutate(sdv_smn1 = ifelse(dpsi_smn1 <= -0.50, T, F)) %>% 
  group_by(category) %>% 
  summarise(num_sdv_smn1 = length(which(sdv_smn1 == T)),
            category_num_smn1 = n()) %>% 
  mutate(percent_sdv_smn1 = num_sdv_smn1/ category_num_smn1) %>%
  arrange(desc(num_sdv_smn1)) %>%
  ungroup()

# percent splice-disrupting mutants by category
snv_by_category_smn1 <- data %>% 
  filter(nat_index_smn1 >= 0.75) %>%
  mutate(snv_smn1 = ifelse(num_changes == 1, T, F)) %>% 
  group_by(category) %>% 
  summarise(num_snv_smn1 = length(which(snv_smn1 == T)),
            category_num_smn1 = n()) %>% 
  # mutate(percent_sdv_smn1 = num_sdv_smn1/ category_num_smn1) %>%
  arrange(desc(num_snv_smn1)) %>%
  ungroup()

sdv_snv_by_category_smn1 <- data %>% 
  filter(nat_index_smn1 >= 0.75, num_changes == 1) %>%
  mutate(sdv_smn1 = ifelse(dpsi_smn1 <= -0.50, T, F)) %>% 
  group_by(category) %>% 
  summarise(num_sdv_smn1 = length(which(sdv_smn1 == T)),
            category_num_smn1 = n()) %>% 
  mutate(percent_sdv_smn1 = num_sdv_smn1/ category_num_smn1) %>%
  arrange(desc(num_sdv_smn1)) %>%
  ungroup()

gg <- data %>% 
  filter(nat_index_smn1 >= 0.75, !is.na(category)) %>%
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
  mutate(percent_sdv_smn1 = num_sdv_smn1/ category_num_smn1 * 100) %>%
  ggplot(aes(x = category_fctr, y = percent_sdv_smn1)) +
  # geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
  #             # aes(color = nat_index_smn1)
  # ) + 
  geom_col(alpha = 1, width = 0.8, fill = "blue") +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 105) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  # breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
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

ggsave(paste0('../../figs/sre/smn1/', 
              'sre_smn1_all_categories_percent_sdv', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = 100)


data %>% 
  filter(nat_index_smn1 >= 0.75, !is.na(category), num_changes == 1) %>%
  dplyr::select(category)

gg1 <- data %>% 
  filter(nat_index_smn1 >= 0.75, !is.na(category), num_changes == 1) %>%
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
  mutate(percent_sdv_smn1 = num_sdv_smn1/ category_num_smn1 * 100) %>%
  ggplot(aes(x = category_fctr, y = percent_sdv_smn1)) +
  # geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
  #             # aes(color = nat_index_smn1)
  # ) + 
  geom_col(alpha = 1, width = 0.8, fill = "blue") +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 105) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  # breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 14, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = "% SDV")

# gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/smn1/', 
              'sre_smn1_all_categories_percent_sdv_single_bp', plot_format), 
       gg1, width = 12, height = 6, dpi = 100)

gg2 <- data %>% 
  filter(nat_index_smn1 >= 0.75, !is.na(category), num_changes == 1) %>%
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
  ggplot(aes(x = category_fctr, y = num_sdv_smn1)) +
  # geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
  #             # aes(color = nat_index_smn1)
  # ) + 
  geom_col(alpha = 1, width = 0.8, fill = "blue") +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 15) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  # breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 14, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = "Number of SDVs")

# gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/smn1/', 
              'sre_smn1_all_categories_num_sdv_single_bp', plot_format), 
       gg2, width = 12, height = 6, dpi = 100)

gg3 <- plot_grid(gg1, gg2, ncol = 1, labels = "AUTO")

ggsave(paste0('../../figs/sre/smn1/', 
              'sre_smn1_sdv_single_bp_both', plot_format), 
       gg3, width = 12, height = 10, dpi = 100)
  # arrange(desc(num_sdv_smn1)) %>%
  # ungroup()

#DHFR SNVs that are SDVs
gg1 <- data %>% 
  filter(nat_index_dhfr >= 0.75, !is.na(category), num_changes == 1) %>%
  mutate(sdv_dhfr = ifelse(dpsi_dhfr <= -0.50, T, F)) %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>%
  group_by(category_fctr) %>% 
  summarise(num_sdv_dhfr = length(which(sdv_dhfr == T)),
            category_num_dhfr = n()) %>% 
  mutate(percent_sdv_dhfr = num_sdv_dhfr/ category_num_dhfr * 100) %>%
  ggplot(aes(x = category_fctr, y = percent_sdv_dhfr)) +
  # geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
  #             # aes(color = nat_index_smn1)
  # ) + 
  geom_col(alpha = 1, width = 0.8, fill = "blue") +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 105) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  # breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 14, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = "% SDV")

# gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/dhfr/', 
              'sre_dhfr_all_categories_percent_sdv_single_bp', plot_format), 
       gg1, width = 12, height = 6, dpi = 100)

gg2 <- data %>% 
  filter(nat_index_dhfr >= 0.75, !is.na(category), num_changes == 1) %>%
  mutate(sdv_dhfr = ifelse(dpsi_dhfr <= -0.50, T, F)) %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>%
  group_by(category_fctr) %>% 
  summarise(num_sdv_dhfr = length(which(sdv_dhfr == T)),
            category_num_dhfr = n()) %>% 
  ggplot(aes(x = category_fctr, y = num_sdv_dhfr)) +
  # geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
  #             # aes(color = nat_index_smn1)
  # ) + 
  geom_col(alpha = 1, width = 0.8, fill = "blue") +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 15) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  # breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 14, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = "Number of SDVs")

# gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/dhfr/', 
              'sre_dhfr_all_categories_num_sdv_single_bp', plot_format), 
       gg2, width = 12, height = 6, dpi = 100)

gg3 <- plot_grid(gg1, gg2, ncol = 1, labels = "AUTO")

ggsave(paste0('../../figs/sre/dhfr/', 
              'sre_dhfr_sdv_single_bp_both', plot_format), 
       gg3, width = 12, height = 10, dpi = 100)
# mutate(category_fctr = factor(category,
#                               levels = c(splice_site_categories, 
#                                          esr_categories, 
#                                          random_exon_categories,
#                                          intron_categories, 
#                                          random_intron_categories, 
#                                          other_categories)))

# DHFR intron backbone INCLUDED
gg <- data %>% 
  filter(seq_type == 'mut', nat_index_dhfr >= 0.5) %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_dhfr))+ 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
              # aes(color = nat_index_smn1)
  ) + 
  geom_boxplot(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = -1), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"), 
        axis.text.y = element_text(size = 10, color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/dhfr/',
              'sre_dhfr_all_categories_no_x_axis_INCLUDED', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = hi_res)

# DHFR intron backbone EXCLUDED
gg <- data %>% 
  filter(seq_type == 'mut', nat_index_dhfr < 0.5) %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50'
              # aes(color = nat_index_smn1)
  ) + 
  geom_boxplot(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = -1), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"), 
        axis.text.y = element_text(size = 10, color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/sre/dhfr/',
              'sre_dhfr_all_categories_no_x_axis_EXCLUDED', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = hi_res)

# percent splice-disrupting mutants by category
sdv_by_category_dhfr <- data %>% 
  filter(nat_index_dhfr >= 0.75) %>%
  mutate(sdv_dhfr = ifelse(dpsi_dhfr <= -0.50, T, F)) %>% 
  group_by(category) %>% 
  summarise(num_sdv_dhfr = length(which(sdv_dhfr == T)),
            category_num_dhfr = n()) %>% 
  mutate(percent_sdv_dhfr = num_sdv_dhfr/ category_num_dhfr) %>%
  arrange(desc(num_sdv_dhfr)) %>%
  ungroup()

# percent splice-disrupting mutants by category
snv_by_category_dhfr <- data %>% 
  filter(nat_index_dhfr >= 0.75) %>%
  mutate(snv_dhfr = ifelse(num_changes == 1, T, F)) %>% 
  group_by(category) %>% 
  summarise(num_snv_dhfr = length(which(snv_dhfr == T)),
            category_num_dhfr = n()) %>% 
  # mutate(percent_sdv_dhfr = num_sdv_dhfr/ category_num_dhfr) %>%
  arrange(desc(num_snv_dhfr)) %>%
  ungroup()

sdv_snv_by_category_dhfr <- data %>% 
  filter(nat_index_dhfr >= 0.75, num_changes == 1) %>%
  mutate(sdv_dhfr = ifelse(dpsi_dhfr <= -0.50, T, F)) %>% 
  group_by(category) %>% 
  summarise(num_sdv_dhfr = length(which(sdv_dhfr == T)),
            category_num_dhfr = n()) %>% 
  mutate(percent_sdv_dhfr = num_sdv_dhfr/ category_num_dhfr) %>%
  arrange(desc(num_sdv_dhfr)) %>%
  ungroup()



# range of number of changes by category
num_changes_summary <- data %>% 
  group_by(category) %>% 
  summarise(min_change = min(num_changes),
            max_change = max(num_changes)) %>%
  ungroup()

###############################################################################
# Exonic motifs
###############################################################################

# calculate change in average HAL score
data <- data %>% 
  dplyr::rename('avg_HAL_score' = 'avg_exon_effect_score') %>% 
  group_by(ensembl_id) %>% 
  mutate(delta_avg_HAL_score = 
           avg_HAL_score - avg_HAL_score[sub_id == '000']) %>% 
  ungroup()

################################
# HAL (hexamer additive linear)
# PMID: 26496609
################################
data <- data %>%
  mutate(HAL_bin = case_when(.$delta_avg_HAL_score < 0 ~ 'down',
                             .$delta_avg_HAL_score > 0 ~ 'up',
                             .$delta_avg_HAL_score == 0 ~ 'same'))

data %>% 
    filter(seq_type == 'mut', HAL_bin != 'same') %>% 
    wilcox.test(dpsi_smn1 ~ HAL_bin, data = .)

# SMN1 intron backbone INCLUDED
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut', nat_index_smn1 >= 0.5) %>%
  ggplot(aes(HAL_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' avg. exon hexamer score')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T, 
                        tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme( axis.title.x = element_text(size = 16, vjust = -2, hjust = 1), 
         axis.title.y = element_text(size = 18, vjust = 40),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_line(color = "grey50"),
         axis.line.x = element_line(color = "grey50"),
         axis.line.y = element_line(color = "grey50"),
         axis.text.x = element_text(size = 18, color = "grey20"), 
         axis.text.y = element_text(color = "grey20"), 
         plot.margin = unit(c(1,1,1,1.5),"mm"),
         legend.position = 'none')

ggsave(paste0('../../figs/sre/smn1/sre_smn1_hal_INCLUDED', 
              plot_format_main), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# SMN1 intron backbone EXCLUDED
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut', nat_index_smn1 < 0.5) %>%
  ggplot(aes(HAL_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' avg. exon hexamer score')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T, 
                        tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme( axis.title.x = element_text(size = 16, vjust = -2, hjust = 1), 
         axis.title.y = element_text(size = 18, vjust = 40),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_line(color = "grey50"),
         axis.line.x = element_line(color = "grey50"),
         axis.line.y = element_line(color = "grey50"),
         axis.text.x = element_text(size = 18, color = "grey20"), 
         axis.text.y = element_text(color = "grey20"), 
         plot.margin = unit(c(1,1,1,1.5),"mm"),
         legend.position = 'none')

ggsave(paste0('../../figs/sre/smn1/sre_smn1_hal_EXCLUDED', 
              plot_format_main), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# DHFR intron backbone
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut', nat_index_dhfr >= 0.5) %>%
  ggplot(aes(HAL_bin, dpsi_dhfr))  + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' avg. exon hexamer score (HAL)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T, 
                        tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 14, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        plot.margin = unit(c(1,1,1,1.5),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/sre/dhfr/sre_dhfr_hal_INCLUDED', 
              plot_format), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# DHFR intron backbone EXCLUDED
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut', nat_index_dhfr < 0.5) %>%
  ggplot(aes(HAL_bin, dpsi_dhfr))  + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' avg. exon hexamer score (HAL)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T, 
                        tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 14, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        plot.margin = unit(c(1,1,1,1.5),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/sre/dhfr/sre_dhfr_hal_EXCLUDED', 
              plot_format), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

########################
# Ke 2011 ESE/ESS 6-mer
# PMID: 21659425
########################
data <- data %>% 
  group_by(ensembl_id) %>% 
  mutate(delta_Ke2011_avg_score = 
           Ke2011_avg_score - Ke2011_avg_score[sub_id == '000'],
         delta_Ke2011_ESE_avg_score = 
           Ke2011_ESE_avg_score - Ke2011_ESE_avg_score[sub_id == '000'],
         delta_Ke2011_ESS_avg_score = 
           Ke2011_ESS_avg_score - Ke2011_ESS_avg_score[sub_id == '000']) %>% 
  ungroup()

data <- data %>%
  mutate(Ke_bin = case_when(.$delta_Ke2011_avg_score < 0 ~ 'down',
                            .$delta_Ke2011_avg_score > 0 ~ 'up',
                            .$delta_Ke2011_avg_score == 0 ~ 'same'),
         Ke_ESE_bin = case_when(.$delta_Ke2011_ESE_avg_score < 0 ~ 'down',
                                .$delta_Ke2011_ESE_avg_score > 0 ~ 'up',
                                .$delta_Ke2011_ESE_avg_score == 0 ~ 'same'),
         Ke_ESS_bin = case_when(.$delta_Ke2011_ESS_avg_score < 0 ~ 'down',
                                .$delta_Ke2011_ESS_avg_score > 0 ~ 'up',
                                .$delta_Ke2011_ESS_avg_score == 0 ~ 'same'))

# SMN1 intron backbone INCLUDED
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut', nat_index_smn1 >= 0.5) %>%
  ggplot(aes(Ke_bin, dpsi_smn1)) +
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' avg. exon hexamer score (Ke)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', 
                        map_signif_level = T, tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 15, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"),
        plot.margin = unit(c(1,1,1,1),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/sre/smn1/sre_smn1_Ke11_INCLUDED', plot_format),
       gg, width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# SMN1 intron backbone EXCLUDED
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut', nat_index_smn1 < 0.5) %>%
  ggplot(aes(Ke_bin, dpsi_smn1)) +
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' avg. exon hexamer score (Ke)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', 
                        map_signif_level = T, tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 15, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"),
        plot.margin = unit(c(1,1,1,1),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/sre/smn1/sre_smn1_Ke11_EXCLUDED', plot_format),
       gg, width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# DHFR intron backbone INCLUDED
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut', nat_index_dhfr >= 0.5) %>%
  ggplot(aes(Ke_bin, dpsi_dhfr))  +
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' avg. exon hexamer score (Ke)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', 
                        map_signif_level = T, tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 15, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"),
        plot.margin = unit(c(1,1,1,1),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/sre/dhfr/sre_dhfr_Ke11_INCLUDED', plot_format),
       gg, width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# DHFR intron backbone EXCLUDED
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut', nat_index_dhfr < 0.5) %>%
  ggplot(aes(Ke_bin, dpsi_dhfr))  +
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' avg. exon hexamer score (Ke)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', 
                        map_signif_level = T, tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 15, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"),
        plot.margin = unit(c(1,1,1,1),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/sre/dhfr/sre_dhfr_Ke11_EXCLUDED', plot_format),
       gg, width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

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

gg <- data %>%
  filter(acc_score_fold_change != 0) %>%
  mutate(acc_fold_change_bin = 
           cut(acc_score_fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(acc_fold_change_bin)) %>%
  ggplot(aes(acc_fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  geom_boxplot(alpha = 0) + 
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = c(0.2, 0.4, 0.6, 0.8, 1), 
                         colors = pal(321)) +
  labs(x = '', y = '',
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = 12, color = "grey20"), 
        text = element_text(size = 12),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.8, "cm"),
        legend.text = element_text(size = 24, color = "grey20"),
        legend.title = element_text(size = 32),
        plot.margin = unit(c(0,0,0,0),"mm"))

legend <- g_legend(gg)
tiff(paste0('../../figs/sre/both/legend', plot_format_main), 
     width = 40, height = 125, units = 'mm', res = hi_res)
grid.newpage()
grid.draw(legend)
dev.off()

legend <- g_legend(gg)
png(paste0('../../figs/sre/both/legend', plot_format), 
     width = 40, height = 125, units = 'mm', res = hi_res)
grid.newpage()
grid.draw(legend)
dev.off()
