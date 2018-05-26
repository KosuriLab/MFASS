load_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
  if(length(new_pkgs)) install.packages(new_pkgs)
  for(pkg in pkgs){
    suppressWarnings(suppressMessages(library(pkg, character.only = T)))
  }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'grid', 'gtable', 'ggsignif', 
          'LaCroixColoR')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.tiff'
hi_res <- 600

data <- read.table('../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T) %>% 
    filter(!is.na(v2_dpsi))

data <- data %>%
    mutate(AF_bin = cut(AF, 
                        breaks = c(0, 0.00001, 0.000025, 0.0005, 0.005, 1), 
                        include.lowest = T,
                        labels = c('Singleton', 'AC = 2-3', 
                                   'AC = 4-10', '0.05-0.5%', '>0.5%' ))) %>% 
    mutate(sdv = ifelse((v2_dpsi <= -0.50) & (nat_v2_index >= 0.5) , 'SDV', 'non-SDV'))

mutant <- data %>% 
  filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

###############################################################################
# Figure 3A, ExAC global allele frequency 
###############################################################################
sdv_mutant_v2 <- mutant %>% 
  filter(nat_v2_index >= 0.50, category == 'mutant') %>%
  group_by(category, AF_bin) %>% 
  summarise(num_sdv = length(which(sdv == 'SDV')),
            category_num = n()) %>% 
  mutate(percent_sdv = num_sdv/ category_num * 100) %>%
  arrange(desc(num_sdv)) %>%
  ungroup()

sdv_mutant_v2 %>% 
  filter(!is.na(AF_bin)) %>% 
  ggplot(aes(AF_bin, percent_sdv)) + 
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

ggsave(paste0("../../figs/snv/fig3A_allele_frequency_percent_sdv", plot_format), 
       width = 4, height = 2.25, units = 'in', dpi = 300)


# reviewer #1: SDV AF distribution vs non-SDV AF
mutant %>% 
    filter(!is.na(AF_bin)) %>% 
    ggplot(aes(AF_bin)) + 
        geom_bar(aes(fill = sdv), position = 'dodge') +
    labs(x = '', fill = 'SDV') +
    scale_fill_manual(values = lacroix_palette('PassionFruit', type='discrete')) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

# SDV vs all
num_all_var <- mutant %>% filter(!is.na(AF_bin)) %>% nrow()
gg1 <- mutant %>% 
    filter(!is.na(AF_bin)) %>% 
    group_by(AF_bin) %>% 
    summarise(percent_AF_bin = (n() / num_all_var) * 100 ) %>% 
    ggplot(aes(AF_bin, percent_AF_bin)) + geom_bar(stat = 'identity') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) + 
    labs(x = '', y = 'percentage', title = 'All ExAC variants') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

num_sdv <- mutant %>% filter(!is.na(AF_bin), sdv == 'SDV') %>% nrow()
gg2 <- mutant %>% 
    filter(!is.na(AF_bin), sdv == 'SDV') %>% 
    group_by(AF_bin) %>% 
    summarise(percent_AF_bin = (n() / num_sdv) * 100 ) %>% 
    ggplot(aes(AF_bin, percent_AF_bin)) + geom_bar(stat = 'identity') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) + 
    labs(x = '', y = 'percentage', title = 'SDV ExAC variants') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

num_non_sdv <- mutant %>% filter(!is.na(AF_bin), sdv == 'non-SDV') %>% nrow()
gg3 <- mutant %>% 
    filter(!is.na(AF_bin), sdv == 'non-SDV') %>% 
    group_by(AF_bin) %>% 
    summarise(percent_AF_bin = (n() / num_non_sdv) * 100 ) %>% 
    ggplot(aes(AF_bin, percent_AF_bin)) + geom_bar(stat = 'identity') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) + 
    labs(x = '', y = 'percentage', title = 'non-SDV ExAC variants') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

plot_grid(gg1, gg2, gg3, nrow = 1)

# SDV vs. non-SDV percentage by functional category
df1 <- mutant %>% 
    filter(!is.na(AF_bin)) %>% 
    filter(sdv == 'SDV') %>%
    group_by(consequence) %>% 
    tally() %>% 
    mutate(pct = n / sum(n),
           sdv = 'SDV')

df2 <- mutant %>% 
    filter(!is.na(AF_bin)) %>% 
    filter(sdv == 'non-SDV') %>%
    group_by(consequence) %>% 
    tally() %>% 
    mutate(pct = n / sum(n),
           sdv = 'non-SDV')

bind_rows(df1, df2) %>% 
    ggplot(aes(consequence, pct*100)) + 
    geom_bar(stat = 'identity', aes(fill = sdv), position = 'dodge') +
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
    scale_fill_manual(values =  c('black', 'red')) +
    labs(x = '', y = 'percentage (within group)', fill = '') +
   coord_flip()

# number of unique genes for mutants
length(unique(mutant$Gene))

# number of unique genes for SDVs
length(unique(filter(mutant, sdv == 'SDV')$Gene))

# per exon proportion of SDVs
mutant %>% 
    select(sdv, exon_id_new) %>% 
    group_by(exon_id_new) %>% 
    summarise(pct = length(which(sdv == 'SDV')) / n()) %>% 
    ggplot(aes(pct*100)) + 
    geom_histogram(binwidth = 1) +
    scale_y_log10() + annotation_logticks(sides = 'l') + 
    labs(x = 'percentage', title = 'Percentage of SDVs per exon')

# gnomad

