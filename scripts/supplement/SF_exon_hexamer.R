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

data %>% 
    filter(seq_type == 'mut', HAL_bin != 'same') %>% 
    wilcox.test(dpsi_smn1 ~ HAL_bin, data = .)

# SMN1 intron backbone INCLUDED
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut', nat_index_smn1 >= 0.5) %>%
  ggplot(aes(HAL_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'grey50') +
  geom_violin(alpha = 0, color = "grey35") +
  # scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), 
  #                        colors = pal(321)) + 
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
  # scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), 
  #                        colors = pal(321)) + 
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
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
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
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
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
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
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
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
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
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
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
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
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
  # scale_colour_gradientn(limits = c(-0.005, 1), 
  #                        breaks = c(0.2, 0.4, 0.6, 0.8, 1), 
  #                        colors = pal(321)) +
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
