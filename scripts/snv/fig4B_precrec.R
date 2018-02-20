### Graph PR and ROC curves for various external models ###

load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'grid', 'pROC')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format_main <- '.tiff'
plot_format <- '.png'
hi_res <- 600

color_cadd = '#FF9912'
color_dann = '#941494'
color_fathmm = '#108C44'
color_fitcons = '#6AA5CD'
color_linsight = '#ABB9B9'
color_spanr =  '#ED1E24'
color_hal = '#000080'
color_phastcons <- 'black'
color_phylop <- 'turquoise'

## ROC curves ##

data <- read.table('../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
spanr <- read.table('../../processed_data/snv/snv_SPANR_scores_capped.txt',
                    sep = '\t', header = T)
hal <- read.table('../../processed_data/snv/snv_HAL_scores.txt',
                  sep = '\t', header = T)

data <- data %>% 
    left_join(select(spanr, id, dpsi_spanr_capped), by = 'id') %>% 
    left_join(select(hal, id, DPSI_pred), by = 'id') %>% 
    filter(category != 'control', nat_v2_index >= 0.5, (!is.na(v2_dpsi))) %>% 
    mutate(sdv = ifelse(strong_lof == T, 1, 0)) %>% 
    filter(!is.na(sdv))

run_roc <- function(df, type) {
    roc_objs <- list()
    if(type == 'exon'){
        roc_objs[['hal']] <- roc(df$sdv, df$DPSI_pred, ci = T)
    }
    roc_objs[['spanr']] <- roc(df$sdv, df$dpsi_spanr_capped, ci = T)
    roc_objs[['dann']] <- roc(df$sdv, df$dann_score, ci = T)
    roc_objs[['linsight']] <- roc(df$sdv, df$linsight_score, ci = T)
    roc_objs[['fitcons']] <- roc(df$sdv, df$fitCons_score, ci = T)
    roc_objs[['fathmm']] <- roc(df$sdv, df$coding_score, ci = T)
    roc_objs[['cadd']] <- roc(df$sdv, df$cadd_score, ci = T)
    roc_objs[['phastCons']] <-roc(data$sdv, data$mean_cons_score, ci = T)
    roc_objs[['phylop']] <- roc(data$sdv, data$phylop_score, ci = T)
    return(roc_objs)
}

all_roc_objs <- run_roc(data, type = 'all')
intron_roc_objs <- run_roc(filter(data, label != 'exon'), type = 'intron')
exon_roc_objs <- run_roc(filter(data, label == 'exon'), type = 'exon')

plot_rocs <- function(roc_objs, exon = F) {
    plot(roc_objs[['spanr']], col = color_spanr, legacy.axes = T, xlab = '', ylab = '')  # changed
    plot(roc_objs[['dann']], add = T, col = color_dann)
    plot(roc_objs[['linsight']], add = T, col = color_linsight)
    plot(roc_objs[['fitcons']], add = T, col = color_fitcons)
    plot(roc_objs[['fathmm']], add = T, col = color_fathmm)
    plot(roc_objs[['cadd']], add = T, col = color_cadd)
    plot(roc_objs[['phastCons']], add = T, col = color_phastcons)
    plot(roc_objs[['phylop']], add = T, col = color_phylop)
    if(exon){
        plot(roc_objs[['hal']], add = T, col = color_hal)
    }
}

methods <- c('HAL', 'SPANR', 'DANN', 'LINSIGHT', 'fitCons', 'FATHMM-MKL', 'CADD')
colors <- c(color_hal, color_spanr, color_dann, color_linsight, color_fitcons,
            color_fathmm, color_cadd)

tiff(paste0(dir,"roc_all.tiff"), width = 3.34, height = 3.44, units = 'in', res = 300)
plot_rocs(all_roc_objs, exon = F)
# legend("bottomright", legend = methods[-1], col = colors[-1], lty = 1, cex = 0.5)
dev.off()

tiff(paste0(dir,"roc_intron.tiff"), width = 3.34, height = 3.44, units = 'in', res = 300)
plot_rocs(intron_roc_objs, exon = F)
# legend("bottomright", legend = methods, col = colors, lty = 1, cex = 0.5)
dev.off()

tiff(paste0(dir,"roc_exon.tiff"), width = 3.34, height = 3.44, units = 'in', res = 300)
plot_rocs(exon_roc_objs, exon = T)
# legend("bottomright", legend = methods, col = colors, lty = 1, cex = 0.5)
dev.off()

# ROC for phyloP and phastCons
roc_phastCons <- roc(data$sdv, data$mean_cons_score)
roc_phylop <- roc(data$sdv, data$phylop_score)
plot(roc_phastCons, col = 'red', legacy.axes = T)
plot(roc_phylop, add = T, col = 'black', legacy.axes = T)
legend('bottomright', legend = c('phastCons', 'phyloP'), col = c('red', 'black'),
       lty = 1, cex = 0.5)

# calculate PR curves
### Graph PR and ROC curves for various external models ###

options(stringsAsFactors = F, warn = -1, warnings = -1)

dir <- '../../figs/snv/'
plot_format_main <- '.png'
plot_format <- '.png'
hi_res <- 600
lo_res <- 300

color_cadd = '#FF9912'
color_dann = '#941494'
color_fathmm = '#108C44'
color_fitcons = '#6AA5CD'
color_linsight = '#ABB9B9'
color_spanr =  '#ED1E24'
color_hal = '#000080'
color_phastcons <- 'black'
color_phylop <- 'turquoise'

###############################################################################
# Precision-recall curves
###############################################################################
pr_curve_all <- 
    read.table('../../processed_data/snv/snv_models_pr_curves_all.txt', 
               sep = '\t', header = T) %>% 
    filter(method != 'hal') %>% # only scores exonic variants
    mutate(type = 'all')

pr_curve_exon <- 
    read.table('../../processed_data/snv/snv_models_pr_curves_exon.txt', 
               sep = '\t', header = T) %>% 
    mutate(type = 'exon')

pr_curve_intron <- 
    read.table('../../processed_data/snv/snv_models_pr_curves_intron.txt', 
               sep = '\t', header = T) %>% 
    mutate(type = 'intron')

pr_curve_info <- bind_rows(pr_curve_all, pr_curve_exon) %>% 
    bind_rows(pr_curve_intron)

# all variants
pr_curve_all %>% 
    filter(method != 'fathmm_noncoding') %>% # pick one FATHMM score method
    mutate(method = factor(method, 
                           labels = c('CADD', 'DANN', 'FATHMM-MKL', 'fitCons', 
                                      'LINSIGHT', 'phastCons',
                                      'phyloP', 'SPANR'))) %>% 
    ggplot(aes(recall, precision)) + 
    geom_line(aes(color = method), size = 1.25) + 
    scale_y_log10(breaks = c(3.8, 10, 100), limits = c(-1,100)) +
    annotation_logticks(sides = 'l') +
    labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
    geom_hline(yintercept = 3.8, linetype = 'dashed', color = 'grey40') +
    scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                  color_fitcons, color_linsight, 
                                  color_phastcons, color_phylop, color_spanr)) +
    theme(
        axis.title.x = element_text(size = 19, vjust = -1.5),
        axis.title.y = element_text(size = 19, vjust = -1.5),
        axis.line.x = element_line(color = 'grey30'),
        axis.line.y = element_line(color = 'grey30'),
        axis.ticks = element_line(color = 'grey30'),
        axis.text = element_text(size = 14, color = 'grey20'),
        legend.justification = 'center',
        legend.direction = 'vertical',
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.height = unit(0.25, "inch"),
        legend.key.width = unit(0.6, "inch"))

ggsave(paste0(dir, 'snv_fig4B_snv_pr_curves_with_legend', 
              plot_format_main), 
       height = 4, width = 7, units = 'in', dpi = hi_res)


pr_curve_all %>% 
    filter(method != 'fathmm_noncoding') %>% # pick one FATHMM score method
    mutate(method = factor(method, 
                           labels = c('CADD', 'DANN', 'FATHMM-MKL', 'fitCons', 
                                      'LINSIGHT', 'phastCons', 'phyloP', 'SPANR'))) %>% 
    ggplot(aes(recall, precision)) + 
    geom_line(aes(color = method), size = 1.25) + 
    scale_y_log10(breaks = c(3.8, 10, 100), limits = c(-1,100)) +
    annotation_logticks(sides = 'l') +
    labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 14, color = 'grey20'),
          axis.line.x = element_line(color = 'grey30'),
          axis.line.y = element_line(color = 'grey30'),
          axis.ticks = element_line(color = 'grey30'),
          legend.position = 'none') +
    geom_hline(yintercept = 3.8, linetype = 'dashed', color = 'grey40') +
    scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                  color_fitcons, color_linsight, 
                                  color_phastcons, color_phylop, color_spanr))

ggsave(paste0(dir, 'snv_fig4B_snv_pr_curves_no_legend', 
              plot_format_main), height = 4, width = 4, units = 'in', dpi = hi_res)

# split by intron/exon
pr_curve_info %>%
    filter(method != 'fathmm_noncoding') %>%# pick one FATHMM score method
    mutate(method = factor(method, labels = c('CADD', 'DANN', 'FATHMM-MKL', 
                                              'fitCons', 'HAL', 'LINSIGHT', 
                                              'phastCons', 'phyloP',
                                              'SPANR')),
           type   = factor(type, levels = c('exon', 'intron', 'all'), 
                           labels = c('Exonic SNV', 'Intronic SNV', 'All SNV'))) %>%
    ggplot(aes(recall, precision)) + geom_line(aes(color = method)) +
    # scale_y_log10() +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 3.8, 10, 100)) +
    annotation_logticks(sides = 'l') +
    facet_grid(~ type) +
    scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                  color_fitcons, color_hal, color_linsight, 
                                  color_phastcons, color_phylop, color_spanr)) +
    labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
    geom_hline(yintercept = 3.8, linetype = 'dashed', color = 'grey40') +
    theme(strip.text = element_text(size = 18.5),
          strip.background = element_rect(fill = "#E8E8E8", color = "white"),
          panel.grid = element_blank(),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(fill = NA, color = 'grey50'),
          axis.title.y = element_text(size = 19, vjust = 2.75),
          axis.title.x = element_text(size = 19, vjust = -2.25),
          axis.text.y = element_text(size = 14, color = 'grey30'),
          axis.text.x = element_text(size = 14, color = 'grey30'),
          axis.line.y = element_line(color = 'grey30'),
          axis.line.x = element_line(color = 'grey30'),
          axis.ticks = element_line(color = 'grey30'),
          legend.justification = 'center',
          legend.direction = 'vertical',
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, 'lines'),
          legend.key.height = unit(0.25, "inch"),
          legend.key.width = unit(0.6, "inch")) +
    guides(colour = guide_legend(override.aes = list(size = 3)))

# graphics.off()
# par("mar")
# par(mar=c(1,1,1,1))
ggsave(paste0("../../figs/snv/fig4B_snv_pr_curves_type", plot_format), 
       height = 4.5, width = 13.75, units = 'in', dpi = lo_res)

# split by intron/exon
pr_curve_info %>%
  filter(method != 'fathmm_noncoding') %>%# pick one FATHMM score method
  mutate(method = factor(method, labels = c('CADD', 'DANN', 'FATHMM-MKL', 
                                            'fitCons', 'HAL', 'LINSIGHT', 
                                            'phastCons', 'phyloP',
                                            'SPANR')),
         type   = factor(type, levels = c('exon', 'intron', 'all'), 
                         labels = c('Exonic SNV', 'Intronic SNV', 'All SNV'))) %>%
  ggplot(aes(recall, precision)) + geom_line(aes(color = method)) +
  # scale_y_log10() +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 3.8, 10, 100)) +
  annotation_logticks(sides = 'l') +
  facet_grid(~ type) +
  scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                color_fitcons, color_hal, color_linsight, 
                                color_phastcons, color_phylop, color_spanr)) +
  labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
  geom_hline(yintercept = 3.8, linetype = 'dashed', color = 'grey40') +
  theme(strip.text = element_text(size = 18.5),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(fill = NA, color = 'grey50'),
        axis.title.y = element_text(size = 19, vjust = 2.75),
        axis.title.x = element_text(size = 19, vjust = -2.25),
        axis.text.y = element_text(size = 14, color = 'grey30'),
        axis.text.x = element_text(size = 14, color = 'grey30'),
        axis.line.y = element_line(color = 'grey30'),
        axis.line.x = element_line(color = 'grey30'),
        axis.ticks = element_line(color = 'grey30'),
        legend.position = 'none') +
  guides(colour = guide_legend(override.aes = list(size = 3)))

# graphics.off()
# par("mar")
# par(mar=c(1,1,1,1))
ggsave(paste0("../../figs/snv/fig4b_snv_pr_curves_type", plot_format), 
       height = 4.5, width = 11.5, units = 'in', dpi = lo_res)
