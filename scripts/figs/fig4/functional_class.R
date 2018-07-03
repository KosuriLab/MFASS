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
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

dpsi_threshold <- -0.50

fig_folder <- '../../../figs/fig4/'

mutant <- data %>% 
    filter(category == 'mutant') %>% 
    mutate(sdv = ifelse(v2_dpsi <= -0.50, 'SDV', 'non-SDV'))

# SDV vs. non-SDV percentage by functional category
mutant <- mutant %>% 
    mutate(consequence_simple = case_when(.$consequence == 'splice_donor_variant' ~ 'splice site',
                                          .$consequence == 'splice_acceptor_variant' ~ 'splice site',
                                          .$consequence == 'synonymous_variant' ~ 'synonymous',
                                          .$consequence == 'missense_variant' ~ 'missense',
                                          .$consequence == 'splice_region_variant' ~ 'splice region',
                                          .$consequence == 'intron_variant' ~ 'intron',
                                          TRUE ~ 'NA'))
mutant %>% 
    group_by(sdv, consequence_simple) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(consequence_simple != 'NA') %>% 
    ggplot(aes(consequence_simple, pct*100)) + 
    geom_bar(stat = 'identity', aes(fill = sdv), position = 'dodge') +
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
    scale_fill_manual(values =  c('black', 'darkred')) +
    labs(x = '', y = 'percentage within\nSDV class', fill = '') +
    coord_flip() +
    theme(legend.position = 'top')

mutant %>% 
    group_by(consequence_simple, sdv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(consequence_simple != 'NA') %>% 
    ggplot(aes(consequence_simple, pct*100)) + 
    geom_bar(stat = 'identity', aes(fill = sdv), position = 'dodge') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    scale_fill_manual(values =  c('black', 'darkred')) +
    labs(x = '', y = 'percentage within\neach functional class', fill = '') +
    coord_flip() +
    theme(legend.position = 'top')

mutant %>% 
    filter(consequence_simple != 'NA') %>% 
    ggplot(aes(v2_dpsi)) + geom_density(aes(color = consequence_simple)) +
    labs(x = expression(paste(Delta, ' inclusion index')), color = 'functional class')

mutant %>% 
    filter(consequence_simple != 'NA') %>% 
    ggplot(aes(consequence_simple, v2_dpsi)) + geom_boxplot(outlier.color = 'grey', outlier.size = 0.5) +
    coord_flip() +
    labs(x = 'functional class', y = expression(paste(Delta, ' inclusion index')))
