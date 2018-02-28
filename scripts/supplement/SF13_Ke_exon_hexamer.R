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
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'red') +
  geom_violin(alpha = 0, color = "grey35") +
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

ggsave(paste0('../../figs/supplement/SF13A_sre_smn1_Ke11_INCLUDED', plot_format),
       gg, width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# SMN1 intron backbone EXCLUDED
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut', nat_index_smn1 < 0.5) %>%
  ggplot(aes(Ke_bin, dpsi_smn1)) +
  geom_jitter(alpha = jitter_alpha, size = 0.5, color = 'blue') +
  geom_violin(alpha = 0, color = "grey35") +
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

ggsave(paste0('../../figs/supplement/SF13B_sre_smn1_Ke11_EXCLUDED', plot_format),
       gg, width = 2.5, height = 3, dpi = hi_res, scale = 1.3)
