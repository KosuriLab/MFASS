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
                   sep = '\t', header = T) %>% 
    filter(!is.na(v2_dpsi))

data <- data %>%
    mutate(AF_bin = cut(AF, 
                        breaks = c(0, 0.00001, 0.000025, 0.0005, 0.005, 1), 
                        include.lowest = T,
                        labels = c('Singleton', 'AC = 2-3', 
                                   'AC = 4-10', '0.05-0.5%', '>0.5%' ))) %>% 
    mutate(sdv = ifelse((v2_dpsi <= -0.50) & (nat_v2_index >= 0.5) , 'SDV', 'non-SDV'))

data <- data %>% 
    mutate(AF_bin_gnomad = cut(AF_gnomad, 
                               breaks = c(0, 0.00001, 0.000025, 0.0005, 0.005, 1), 
                               include.lowest = T,
                               labels = c('Singleton', 'AC = 2-3', 
                                          'AC = 4-10', '0.05-0.5%', '>0.5%' )))

mutant <- data %>% 
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

fig_folder <- '../../../figs/fig5/'

###############################################################################
# Figure 5A, gnomAD global allele frequency 
###############################################################################
sdv_mutant_v2 <- mutant %>% 
    filter(nat_v2_index >= 0.50, category == 'mutant') %>%
    group_by(category, AF_bin_gnomad) %>% 
    summarise(num_sdv = length(which(sdv == 'SDV')),
              category_num = n()) %>% 
    mutate(percent_sdv = num_sdv/ category_num * 100) %>%
    arrange(desc(num_sdv)) %>%
    ungroup()

sdv_mutant_v2 %>% 
    filter(!is.na(AF_bin_gnomad)) %>% 
    ggplot(aes(AF_bin_gnomad, percent_sdv)) + 
    geom_col(width = 0.8, color = "#000080", fill = "#000080") + 
    geom_hline(yintercept = 3.8, linetype = "dashed", color = "grey40") +
    labs(x = "ExAC global allele frequency", 
         y = "% SDV") +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = 4.5) +
    theme_bw() + 
    theme(legend.position = 'none', 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(size = 17, vjust = 12),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "grey20"),
          axis.text.x = element_blank()
    ) 

ggsave(paste0(fig_folder, "allele_frequency_gnomad_percent_sdv", plot_format), 
       width = 4, height = 2.25, units = 'in', dpi = 300)


# reviewer #1: SDV AF distribution vs non-SDV AF
mutant %>% 
    filter(!is.na(AF_bin_gnomad)) %>% 
    group_by(sdv, AF_bin_gnomad) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    ggplot(aes(AF_bin_gnomad, pct * 100)) +
    geom_bar(stat = 'identity', position = 'dodge', aes(fill = sdv)) +
    scale_fill_manual(values = c('black', 'darkred')) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
    labs(x = '', y = 'percentage', fill = '') + 
    theme(legend.position = 'top', axis.text.x = element_blank())
ggsave(paste0("../../figs/snv/allele_frequency_sdv_nonsdv", plot_format), 
       width = 5, height = 3.5, units = 'in', dpi = 300)

# SDVs enriched for singletons?
df <- mutant %>% 
    filter(!is.na(AF_bin_gnomad)) %>% 
    mutate(singleton = ifelse(AF_bin_gnomad == 'Singleton', T, F)) %>% 
    group_by(sdv, singleton) %>% 
    tally() %>% 
    ungroup() %>% 
    arrange(desc(sdv), desc(singleton))

mat <- matrix(df$n, nrow = 2, dimnames = list(c('singleton', 'not singleton'),
                                              c('SDV', 'non-SDV')))
fisher.test(mat)
# odds ratio 1.3, more SDV singletons than expected by chance


# # PTV rate
# data <- data %>% 
#     mutate(sdv = ifelse(v2_dpsi <= -0.50, 'SDV', 'non-SDV'),
#            ptv = ifelse(consequence %in% c('stop_gained', 
#                                            'splice_acceptor_variant',
#                                            'splice_donor_variant'),
#                         'PTV', 'non-PTV'))
# counts <- data %>% 
#     filter(category == 'mutant', !is.na(AF_bin)) %>% 
#     group_by(sdv, AF_bin) %>% 
#     tally() %>% 
#     mutate(pct = n / sum(n)) %>% 
#     ungroup() %>% 
#     filter(sdv == 'SDV') %>% 
#     select(-sdv) %>% 
#     mutate(type = 'SDV') 
# 
# counts <- data %>% 
#     filter(category == 'mutant', !is.na(AF_bin)) %>% 
#     group_by(ptv, AF_bin) %>% 
#     tally() %>% 
#     mutate(pct = n / sum(n)) %>% 
#     ungroup() %>% 
#     filter(ptv == 'PTV') %>% 
#     mutate(type = 'PTV') %>% 
#     select(-ptv) %>% 
#     bind_rows(counts)
# 
# ggplot(counts, aes(AF_bin, pct*100)) +
#     geom_bar(stat = 'identity', position = 'dodge', aes(fill = type)) +
#     labs(x = 'allele frequency', y = 'percentage', fill = '') +
#     scale_fill_manual(values = c('darkblue', 'darkred')) +
#     theme(legend.position = 'top')
