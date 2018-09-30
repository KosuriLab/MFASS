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

data <- read.table('../../../processed_data/snv/snv_func_annot.txt', sep = '\t', header = T)
data <- data %>%
    filter(category == 'mutant', nat_v2_index >= 0.5, !(is.na(v2_dpsi)))

fig_folder <- '../../../figs/fig5/'


##############################
# Figure 3b: pLI
##############################

dpsi_threshold <- -0.50
data <- data %>% 
    mutate(sdv = ifelse(v2_dpsi <= dpsi_threshold & nat_v2_index >= 0.50, 'SDV', 'non-SDV'))

counts <- data %>% 
    filter(!is.na(pLI_high), !is.na(sdv), category == 'mutant') %>% 
    mutate(tolerance = ifelse(pLI_high == T, 'intolerant', 'tolerant')) %>% 
    group_by(tolerance, sdv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n))

color_pli = "black"
counts %>%
    filter(sdv == 'SDV') %>%
    ggplot(aes(tolerance, pct*100)) +
    geom_col(width = 0.65, color = color_pli, fill = color_pli) +
    ylab("% SDV") + xlab("pLI") + ylim(0,5) +
    geom_hline(yintercept = 1050/27733*100, linetype = "dashed", color = "grey40") +
    scale_y_continuous(expand = c(0, 0)) +
    expand_limits(y = 4.8) +
    # coord_equal(1/1.125) +
    theme_bw() +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(size = 16, color = "black", vjust = 2),
          axis.title.x = element_text(size = 20, color = "black", vjust = -0.5),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "grey20"),
          axis.text.x = element_text(size = 13, color = "black"))

ggsave(paste0(fig_folder, "pLI_enrichment", plot_format), 
       width = 2.25, height = 3.5, units = 'in', dpi = hi_res)

# current matrix arrangement tests for significantly greater proportion of non-SDVs
# in intolerant genes, equivalent to significantly fewer proportion of SDVS in intolerant
# genes if you rearranged the matrix
df <- matrix(counts$n, nrow = 2)
fisher.test(df, alternative = 'two.sided')


# compare rate of PTVs in pLI genes to SDVs
data <- data %>% 
    mutate(ptv = ifelse(consequence %in% c('stop_gained', 
                                           'splice_acceptor_variant',
                                           'splice_donor_variant'),
                        'PTV', 'non-PTV'))

counts <- data %>% 
    filter(!is.na(pLI_high)) %>% 
    mutate(tolerance = ifelse(pLI_high == T, 'intolerant', 'tolerant')) %>% 
    group_by(tolerance, ptv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    rename(type = ptv) %>% 
    bind_rows(rename(counts, type = sdv))

counts %>% 
    filter(type == 'PTV' | type == 'SDV') %>% 
    ggplot(aes(tolerance, pct*100)) + 
    geom_bar(stat = 'identity', position = 'dodge', aes(fill = type)) + 
    scale_fill_manual(values = c('darkred', 'darkblue')) + 
    labs(x = 'pLI', y = 'percentage', fill = '') + ylim(0,5) +
    geom_hline(yintercept = 1050/27733*100, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 340/27733*100, linetype = "dashed", color = "grey40") +
    scale_y_continuous(expand = c(0, 0)) + 
    expand_limits(y = 4.8) +
    # coord_equal(1/1.125) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(size = 16, color = "black", vjust = 2),
          axis.title.x = element_text(size = 20, color = "black", vjust = -0.5),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "grey20"),
          axis.text.x = element_text(size = 13, color = "black"), 
          legend.position = 'none') 

ggsave(paste0(fig_folder, "pLI_sdv_vs_ptv", plot_format), 
       width = 2.75, height = 4, units = 'in', dpi = hi_res)

fisher.test(matrix(filter(counts, grepl('PTV', type))$n, nrow=2),
            alternative = 'two.sided')
