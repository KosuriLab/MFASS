load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'forcats', 'gridExtra', 
          'grid', 'Unicode', 'svglite')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

hi_res <- 600
plot_format <- '.pdf'

data <- read.table('../../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
    filter((!is.na(v2_dpsi)))

dpsi_threshold <- -0.50

fig_folder <- '../../../figs/fig4/'

data <- data %>% 
    mutate(sdv = ifelse(v2_dpsi <= -0.50 & nat_v2_index >= 0.50, 'SDV', 'non-SDV'))

# SDV vs. non-SDV percentage by functional category
data <- data %>% 
    mutate(consequence_simple = case_when(.$consequence == 'splice_donor_variant' ~ 'splice site',
                                          .$consequence == 'splice_acceptor_variant' ~ 'splice site',
                                          .$consequence == 'synonymous_variant' ~ 'synonymous',
                                          .$consequence == 'missense_variant' ~ 'missense',
                                          .$consequence == 'splice_region_variant' ~ 'splice region',
                                          .$consequence == 'intron_variant' ~ 'intron',
                                          TRUE ~ 'NA'))

# mutant <- filter(mutant, consequence_simple != 'NA')
data$consequence_simple <- factor(data$consequence_simple)
data$consequence_simple <- factor(data$consequence_simple,
                                    levels = levels(data$consequence_simple)[c(2, 6, 4, 5, 1, 3)])

mutant <- data %>% filter(category == 'mutant')
order <- rev(levels(mutant$consequence_simple))[-1]

mutant %>% 
    group_by(consequence_simple, sdv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(consequence_simple != 'NA') %>% 
    ggplot(aes(consequence_simple, pct*100)) + 
    geom_bar(stat = 'identity', aes(fill = sdv), position = 'dodge') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    scale_x_discrete(limits = order) +
    scale_fill_manual(values =  c('black', 'darkblue')) +
    labs(x = '', y = '% within\neach functional class', fill = '') +
    coord_flip() +
    theme(legend.position = 'right')

ggsave(paste0(fig_folder, 'sdv_nonsdv_category_percentage', plot_format),
       width = 5.5, height = 4)

mutant %>% 
    group_by(sdv, consequence_simple) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    ggplot(aes(consequence_simple, pct*100)) + 
    geom_bar(stat = 'identity', aes(fill = sdv), position = 'dodge') +
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
    scale_x_discrete(limits = order) +
    scale_fill_manual(values =  c('black', 'darkblue')) +
    labs(x = '', y = '% overall', fill = '') +
    coord_flip() +
    theme(legend.position = 'none')

ggsave(paste0(fig_folder, 'sdv_nonsdv_overall_percentage', plot_format),
       width = 5, height = 4)


mutant %>% 
    group_by(consequence_simple, sdv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(consequence_simple != 'NA', sdv == 'SDV') %>% 
    ggplot(aes(consequence_simple, pct*100)) + 
    geom_bar(stat = 'identity', position = 'dodge', width = 0.6) +
    scale_x_discrete(limits = order, position = 'top') +
    scale_fill_manual(values =  c('black', 'darkblue')) +
    geom_hline(yintercept = (1050/27733) * 100, linetype = 'dashed', color = 'grey10') +
    labs(x = '', y = 'SDV % within\neach functional class', fill = '') +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 30)) +
    coord_flip() + scale_y_reverse(expand = c(0,0)) +
    theme(legend.position = 'right')

ggsave(paste0(fig_folder, 'sdv_category_percentage', plot_format),
       width = 4.5, height = 5, units = 'in')

mutant %>% 
    group_by(sdv, consequence_simple) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(consequence_simple != 'NA', sdv == 'SDV') %>% 
    ggplot(aes(consequence_simple, pct*100)) + 
    geom_bar(stat = 'identity', position = 'dodge', width = 0.6) +
    scale_x_discrete(limits = order) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 30)) +
    scale_fill_manual(values =  c('black', 'darkblue')) +
    labs(x = '', y = 'SDV % within\neach functional class', fill = '') +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 30),
          plot.margin=unit(c(0.5,1,0,0),"cm")) +
    coord_flip() + 
    theme(legend.position = 'right')

ggsave(paste0(fig_folder, 'sdv_overall_percentage', plot_format),
       width = 4.5, height = 5)

# mutant %>% 
#     filter(consequence_simple != 'NA') %>% 
#     ggplot(aes(v2_dpsi)) + geom_density(aes(color = consequence_simple)) +
#     labs(x = expression(paste(Delta, ' inclusion index')), color = 'functional class')

mutant %>% 
    filter(consequence_simple != 'NA', sdv == 'SDV') %>% 
    ggplot(aes(consequence_simple, v2_dpsi)) + geom_boxplot(outlier.color = 'grey', outlier.size = 0.5) +
    coord_flip() +
    labs(x = 'functional class', y = expression(paste(Delta, ' inclusion index'))) +
    scale_x_discrete(limits = order)

ggsave(paste0(fig_folder, 'sdv_category_boxplot', plot_format),
       width = 5, height = 3.5)
