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

# calculate change in donor/acceptor score between mutant and wild-type 
data <- data %>% 
  group_by(ensembl_id) %>% 
  filter(any(sub_id == '000')) %>% 
  mutate(correct_acc_score_nat = correct_acc_score[sub_id == '000'],
         correct_don_score_nat = correct_don_score[sub_id == '000'],
         delta_acc_score = correct_acc_score - correct_acc_score_nat,
         delta_don_score = correct_don_score - correct_don_score_nat,
         # calculate fold change in score relative to wild-type
         don_score_fold_change = (correct_don_score - correct_don_score_nat) / 
           abs(correct_don_score_nat),
         acc_score_fold_change = (correct_acc_score - correct_acc_score_nat) / 
           abs(correct_acc_score_nat)) %>%
  ungroup()

###############################################################################
# MaxEnt 
###############################################################################

# MaxEnt: splice acceptor and donor score fold-change
# SMN1 intron backbone
# cor.test(data$delta_acc_score, data$dpsi_smn1)
gg <- data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut', nat_index_smn1 >= 0.5) %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  ggplot(aes(fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'red') +
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        axis.title.x = element_text(size = 16, vjust = -2), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 12, color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        text = element_text(size = 12),
        plot.margin = unit(c(0,0,3,0),"mm"),
        legend.position = 'none') 

ggsave(paste0('../../figs/sre/smn1/ext_fig4A_smn1_don_acc_INCLUDED', 
              plot_format), gg, 
       width = 4.6, height = 3, dpi = hi_res, scale = 1.3)

# MaxEnt: splice acceptor and donor score fold-change
# SMN1 intron backbone
# cor.test(data$delta_acc_score, data$dpsi_smn1)
gg <- data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut', nat_index_smn1 < 0.5) %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  ggplot(aes(fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'blue') +
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        axis.title.x = element_text(size = 16, vjust = -2), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 12, color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        text = element_text(size = 12),
        plot.margin = unit(c(0,0,3,0),"mm"),
        legend.position = 'none') 

ggsave(paste0('../../figs/sre/smn1/ext_fig4B_smn1_don_acc_EXCLUDED', 
              plot_format), gg, 
       width = 4.6, height = 3, dpi = hi_res, scale = 1.3)

###############################################################################
# Exonic motifs
###############################################################################

# calculate change in average HAL score
data <- data %>% 
    dplyr::rename(avg_HAL_score = avg_exon_effect_score) %>% 
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



# SMN1 intron backbone INCLUDED
gg <- data %>%
    filter(HAL_bin != 'same', seq_type == 'mut', nat_index_smn1 >= 0.5) %>%
    ggplot(aes(HAL_bin, dpsi_smn1)) + 
    geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'red') +
    geom_violin(alpha = 0, color = "grey35") +
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

ggsave(paste0('../../figs/sre/smn1/ext_fig4C_smn1_hal_INCLUDED', 
              plot_format), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# SMN1 intron backbone EXCLUDED
gg <- data %>%
    filter(HAL_bin != 'same', seq_type == 'mut', nat_index_smn1 < 0.5) %>%
    ggplot(aes(HAL_bin, dpsi_smn1)) + 
    geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'blue') +
    geom_violin(alpha = 0, color = "grey35") +
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

ggsave(paste0('../../figs/sre/smn1/ext_fig4D_smn1_hal_EXCLUDED', 
              plot_format), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)
