load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'grid', 'gtable', 'ggsignif')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.pdf'
hi_res <- 600

data <- read.table('../../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

fig_folder <- '../../../figs/fig5/'

############################################################
# Fig. 5c, phyloP: SDVs vs non-SDVs
############################################################
#phyloP summary
phylopSummary <- data %>%
    dplyr::group_by(strong_lof) %>%
    dplyr::summarize(phylop_median = median(phylop_score),
                     phylop_se = sqrt(var(phylop_score)/length(phylop_score)))

data %>%
    ggplot() +
    geom_violin(aes(strong_lof, phylop_score, fill = strong_lof, color = strong_lof)) +
    geom_boxplot(aes(strong_lof, phylop_score), width = 0.05, outlier.shape=NA) +
    theme_bw() + 
    theme(legend.position = 'none',
          strip.text = element_text(size = 18.5),
          strip.background = element_rect(fill = "#E0E0E0", color = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          plot.margin = unit(c(0,0,0,0),'in')) +
    scale_x_discrete(labels = c('non-\nSDV', 'SDV')) +
    ylab('phyloP score')

ggsave(filename = paste0(fig_folder, "snv_phylop_violin", plot_format),
       width = 2, height = 3, units = 'in') 

# could missense variants be driving the difference in phyloP scores?
with_missense <- data %>% filter(!is.na(strong_lof))
wilcox.test(phylop_score ~ strong_lof, with_missense)
without_missense <- data %>% filter(!is.na(strong_lof), consequence != 'missense_variant')
wilcox.test(phylop_score ~ strong_lof, without_missense)


############################################################
# Fig. 3DE, phyloP: proportion of SDVs across gene regions
############################################################

data <- data %>% 
    mutate(cons_bin = cut(phylop_score, breaks = c(-10, -2.0, -1.2, 1.2, 2.0, 10), 
                          include.lowest = T, 
                          labels = c('accelerating', 'drop1', 'neutral', 
                                     'drop2', 'deleterious')))
data$label_renamed <- factor(data$label, 
                             levels=c("upstr_intron", "exon", "downstr_intron"), 
                             labels=c("Intron\nupstr.", "Exon", "Intron\ndownstr."))    

counts <- data %>% 
    filter(!is.na(strong_lof), !is.na(cons_bin)) %>% 
    filter(cons_bin != 'drop1', cons_bin != 'drop2') %>% 
    group_by(label_renamed, cons_bin, strong_lof) %>% 
    tally() %>% 
    mutate(percent = (n / sum(n)) * 100)
    

counts$cons_bin <- factor(counts$cons_bin, 
                          levels=c('accelerating', 'neutral', 'deleterious'))

counts %>% 
    filter(strong_lof == T) %>% 
    ggplot(aes(label_renamed, percent, fill = factor(cons_bin))) +
    geom_histogram(stat = 'identity', width = 0.75, position = "dodge") +
    ylab("% SDV") +
    xlab("") +
    # facet_wrap(~ cons_bin) +
    geom_hline(yintercept = 3.8, linetype = "dashed", color = "grey20") + 
    scale_y_continuous(breaks = c(0, 3.8, 10, 20, 30, 40), expand = c(0,0)) +
    expand_limits(y = 42.5) +
    theme_bw() + 
    theme(strip.text = element_text(size = 18.5),
          strip.background = element_rect(fill = "#E0E0E0", color = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(size = 19, vjust = 15),
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.4,0.85),
          legend.text = element_text(size = 16),
          legend.key.height = unit(0.25, "inch"),
          legend.key.width = unit(0.4, "inch")
    ) +
    scale_fill_manual(values = c('#458B00','#2E9FFE', 'purple'))

ggsave(paste0(fig_folder, "phyloP_comparison_prop", plot_format),
       width = 4.5, height = 5, units = 'in', dpi = 300)


###############################################################################
# Fig. 5e: absolute counts for LoF SNPs, by phyloP conservation
###############################################################################

counts %>% 
    filter(strong_lof == T) %>% 
    ggplot(aes(label_renamed, n, fill = factor(cons_bin))) +
    geom_histogram(stat = 'identity', width = 0.75, position = "dodge") +
    ylab("Number of SDVs") +
    xlab("") +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = 450) +
    theme_bw() + 
    theme(strip.text = element_text(size = 19),
          strip.background = element_rect(fill = "#E0E0E0", color = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey50"),
          axis.title.y = element_text(size = 19, vjust = 20),
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.ticks.x = element_blank(),
          legend.position = c(0.4,0.85),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.height = unit(0.25, "inch"),
          legend.key.width = unit(0.4, "inch")
    ) +
    scale_fill_manual(values = c('#458B00',
                                 '#2E9FFE',
                                 'purple'))

ggsave(paste0(fig_folder, "phylop_abs_num", plot_format),
       width = 4.5, height = 5, units = 'in', dpi = 300)
