readRenviron("~/.Renviron")
setwd(Sys.getenv("SPLICEVAR"))

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

# calculate PR curves
### Graph PR and ROC curves for various external models ###

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.png'
hi_res <- 300

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
    read.table('../../processed_data/exac/exac_models_pr_curves_all.txt', 
               sep = '\t', header = T) %>% 
    filter(method != 'hal') %>% # only scores exonic variants
    mutate(type = 'all')

pr_curve_exon <- 
    read.table('../../processed_data/exac/exac_models_pr_curves_exon.txt', 
               sep = '\t', header = T) %>% 
    mutate(type = 'exon')

pr_curve_intron <- 
    read.table('../../processed_data/exac/exac_models_pr_curves_intron.txt', 
               sep = '\t', header = T) %>% 
    mutate(type = 'intron')

pr_curve_info <- bind_rows(pr_curve_all, pr_curve_exon) %>% 
    bind_rows(pr_curve_intron)

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
ggsave(paste0("../../figs/exac/exac_suppfig_exac_pr_curves_type", plot_format), 
       height = 4.5, width = 11.5, units = 'in', dpi = hi_res)
