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

plot_format <- '.png'
hi_res <- 600

data <- read.table('../../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

color_unknown = '#FF6103' 
color_probably_damaging = '#002F6C'
color_possibly_damaging = '#6AA5CD'
color_benign = '#9D1309'
color_noannot = '#AAAAAA'

data <- data %>% 
    separate(SIFT, into = c('SIFT', 'SIFT_other'), sep = '[(]', remove = T) %>% 
    separate(PolyPhen, into = c('PolyPhen', 'PolyPhen_other'), sep = '[(]', remove = T) %>% 
    select(-SIFT_other, -PolyPhen_other)

fig_folder <- '../../../figs/fig6/'

data <- data %>% 
    mutate(sdv = ifelse(v2_dpsi <= -0.5, 'SDV', 'non-SDV'))

################
#### Polyphen
################

data %>%
    filter(consequence == "missense_variant", v2_dpsi <= -0.5, !is.na(PolyPhen)) %>%
    ggplot(aes(x = PolyPhen, fill = PolyPhen)) + 
    geom_histogram(stat="count", width = 0.6) + 
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = 180) +
    labs(x = '', #PolyPhen prediction\n for missense SDVs
         y = 'count') +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(vjust = 20, size = 20),
          axis.title.x = element_text(vjust = -20, size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
    ) +
    scale_fill_manual(values = c(color_unknown, 
                                 color_probably_damaging, 
                                 color_possibly_damaging, 
                                 color_benign,
                                 color_noannot)) +
    scale_x_discrete(labels = c('no annotation', 
                                'benign', 
                                'possibly\ndamaging', 
                                'probably\ndamaging', 
                                'unknown')) +
    coord_flip()

ggsave(paste0(fig_folder, 'polyphen_SDVs_missense', plot_format), 
       width = 5, height = 3, dpi = hi_res) 

data %>%
    filter(consequence == "missense_variant", !is.na(PolyPhen)) %>%
    group_by(sdv, PolyPhen) %>% 
    tally() %>% 
    mutate(pct = (n / sum(n)) * 100) %>% 
    ggplot(aes(PolyPhen, pct)) +
    geom_bar(stat = 'identity', aes(fill = sdv), width = 0.6, position = 'dodge') +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = '', 
         y = 'percentage', fill = '') +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey50"),
          # axis.title.y = element_text(vjust = 20, size = 20),
          # axis.title.x = element_text(vjust = -20, size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
    ) +
    scale_fill_manual(values = c('black', 'darkblue')) +
    scale_x_discrete(labels = c('no annotation', 
                                'benign', 
                                'possibly\ndamaging', 
                                'probably\ndamaging', 
                                'unknown')) +
    coord_flip()

ggsave(paste0(fig_folder, 'polyPhen_all_missense', plot_format),
       width = 5, height = 3)

############
#### SIFT
############

color_tolerated =  '#FF6103'
color_deleterious = '#9D1309'
color_noannot = '#AAAAAA'

data %>%
    filter(consequence == "missense_variant", v2_dpsi <= -0.5, !is.na(SIFT)) %>%
    filter(SIFT != "tolerated_low_confidence", SIFT != "deleterious_low_confidence") %>%
    ggplot(aes(x = SIFT, fill = SIFT)) +
    geom_histogram(stat="count", width = 0.6) + 
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = 180) +
    labs(x = '', 
         y = 'count') +
    scale_x_discrete(labels = c('no annotation', 
                                'deleterious', 
                                'tolerated')) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 20, vjust = 60),
          axis.title.x = element_text(vjust = -20, size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14)
    ) +
    coord_flip() +
    scale_fill_manual(values = c(color_tolerated, 
                                 color_deleterious, 
                                 color_noannot))

ggsave(paste0(fig_folder, 'sift_SDVs_missense', plot_format), 
       width = 5, height = 2.1, dpi = hi_res) 


data %>%
    filter(consequence == "missense_variant", !is.na(SIFT)) %>%
    filter(SIFT != "tolerated_low_confidence", SIFT != "deleterious_low_confidence") %>%
    group_by(sdv, SIFT) %>% 
    tally() %>% 
    mutate(pct = (n / sum(n)) * 100) %>% 
    ggplot(aes(SIFT, pct)) +
    geom_bar(stat = 'identity', aes(fill = sdv), width = 0.6, position = 'dodge') +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = '', 
         y = 'percentage', fill = '') +
    scale_fill_manual(values = c('black', 'darkblue')) +
    scale_x_discrete(labels = c('no annotation', 
                                'deleterious', 
                                'tolerated')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.title.y = element_text(size = 20, vjust = 60),
          # axis.title.x = element_text(vjust = -20, size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14)
    ) +
    coord_flip()

ggsave(paste0(fig_folder, 'SIFT_all_missense', plot_format),
       width = 5.5, height = 3)
