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
plot_format <- '.png'
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
cor.test(data$delta_acc_score, data$dpsi_smn1)
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
  # dplyr::select(acc_score_fold_change, don_score_fold_change, splice_site, 
  # fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), 
  #                        colors = pal(321)) +
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

ggsave(paste0('../../figs/sre/smn1/sre_smn1_don_acc_INCLUDED', 
              plot_format_main), gg, 
       width = 4.6, height = 3, dpi = hi_res, scale = 1.3)

# MaxEnt: splice acceptor and donor score fold-change
# SMN1 intron backbone
cor.test(data$delta_acc_score, data$dpsi_smn1)
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
  # dplyr::select(acc_score_fold_change, don_score_fold_change, splice_site, 
  # fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), 
  #                        colors = pal(321)) +
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

ggsave(paste0('../../figs/sre/smn1/sre_smn1_don_acc_EXCLUDED', 
              plot_format_main), gg, 
       width = 4.6, height = 3, dpi = hi_res, scale = 1.3)

# MaxEnt: splice acceptor and donor score fold-change)
# DHFR intron backbone, included

gg <- data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut', nat_index_dhfr >= 0.5) %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  # dplyr::select(acc_score_fold_change, don_score_fold_change, splice_site,
  # fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), 
  #                        colors = pal(321)
  #                        ) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        axis.title.x = element_text(size = 16, vjust = -2), 
        axis.title.y = element_text(size = 18), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text = element_text(size = 12, color = "grey20"), 
        text = element_text(size = 12),
        plot.margin = unit(c(0,0,3,0),"mm"),
        legend.position = 'none') 

ggsave(paste0('../../figs/sre/dhfr/sre_dhfr_don_acc_INCLUDED', 
              plot_format_main), 
       gg, width = 4.6, height = 3, dpi = hi_res, scale = 1.3)

# MaxEnt: splice acceptor and donor score fold-change)
# DHFR intron backbone, included

gg <- data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut', nat_index_dhfr < 0.5) %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  # dplyr::select(acc_score_fold_change, don_score_fold_change, splice_site,
  # fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), 
  #                        colors = pal(321)) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        axis.title.x = element_text(size = 16, vjust = -2), 
        axis.title.y = element_text(size = 18), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text = element_text(size = 12, color = "grey20"), 
        text = element_text(size = 12),
        plot.margin = unit(c(0,0,3,0),"mm"),
        legend.position = 'none') 

ggsave(paste0('../../figs/sre/dhfr/sre_dhfr_don_acc_EXCLUDED', 
              plot_format_main), 
       gg, width = 4.6, height = 3, dpi = hi_res, scale = 1.3)

