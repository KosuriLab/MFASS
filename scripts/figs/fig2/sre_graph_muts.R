library(dplyr)
library(ggplot2)
library(cowplot)

options(stringsAsFactors = F, scipen = 10000)

data <- read.table('../../../processed_data/sre/sre_data_mut_changes.txt', sep = '\t', header = T)

data <- data %>% 
    mutate(intron_changes = intron1_changes + intron2_changes)

data$category_rename <- factor(data$category)
levels(data$category_rename) <- c('all_mut_both', 'all_mut_exon', 'all_mut_intron',
                                  'strongest_ESE', 'strongest_ESS', 'strongest_AICS',
                                  'strongest_DICS', 'consrv_1nt', 'consrv_3nt',
                                  'weaken_spurious_acceptor', 'weaken_spurious_donor', 'destroy_acceptor',
                                  'destroy_donor', 'all_mut_random_intron', 'weaken_acceptor_donor_region',
                                  'rbp_site', 'all_ESE', 'all_ESS',
                                  'destroy_spurious_acceptor', 'destroy_spurious_donor',
                                  'all_AICS', 'all_DICS', 'rnd_exon_1nt', 'rnd_exon_2nt',
                                  'rnd_exon_3nt', 'rnd_exon_5nt', 'rnd_intron_1nt',
                                  'rnd_intron_2nt', 'rnd_intron_3nt', 'rnd_intron_5nt',
                                  'same_score_acceptor', 'same_score_donor', 'variation',
                                  'weaken_acceptor_region', 'weaken_donor_region')

data <- data %>% 
    mutate(sdv = ifelse(dpsi_smn1 <= -0.50 & nat_index_smn1 >= 0.50, 'SDV', 'non-SDV'))


data$category_rename <- factor(data$category_rename,
                                       levels = levels(data$category_rename)[c(34, 35, 15, 12, 13, 31, 32, 10, 11, 19, 20,
                                                               4, 17, 5, 18, 23:26, 2, 
                                                               6, 7, 21, 22, 27:30, 3, 14,
                                                               1, 8, 9, 16, 33)])

gg1 <- data %>% 
    filter(!is.na(category_rename), seq_type == 'mut') %>% 
    ggplot(aes(category_rename, dpsi_smn1)) + 
        geom_boxplot(outlier.fill = 'grey', outlier.size = 0.5, outlier.color = 'grey') + 
        coord_flip() +
        geom_hline(yintercept = 0, linetype = 'dashed', color= 'grey23') +
        geom_hline(yintercept = -0.5, linetype = 'dashed', color= 'grey23') +
        scale_x_discrete(limits = rev(levels(data$category_rename)),
                         position = 'top') +
    labs(x = '', y = expression(paste(Delta, ' inclusion index'))) +
    scale_y_continuous(breaks = seq(-1, 1, 0.50)) +
    theme(axis.text.y = element_blank())
gg1

ggsave(filename = '../../figs/sre/category_delta_index.pdf', width = 3.5, height = 7)

num_changes <- data %>% 
    filter(seq_type == 'mut', !is.na(sdv)) %>% 
    group_by(category_rename) %>% 
    summarise(mean_change = mean(num_changes),
              sd_change = sd(num_changes))

percent_sdv <- data %>%  
    filter(seq_type == 'mut', !is.na(sdv)) %>% 
    group_by(category_rename, sdv) %>% 
    tally() %>% 
    mutate(pct = 100 * (n / sum(n))) %>% 
    filter(sdv == 'SDV')
    


