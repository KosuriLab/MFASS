### Process raw data with appropriate normalization and filters ###
load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs, repos = 'http://cran.stat.ucla.edu/')
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'weights', 'wCorr')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.tiff'

###############################################################################
# Read in data
###############################################################################

# snv first sequencing run (v1), three bins
# select sort 2 only
snv_v1 <- read.csv('../../processed_data/snv/snv_v1_all_alignments.csv') %>% 
    select(-ends_with('S1'), id,
           DP_R1 = DP1S2, DP_R2 = DP2S2, INT_R1 = INT1S2, INT_R2 = INT2S2,
           PS_R1 = R1.PS, PS_R2 = R2.PS, SP_R1 = SP1S2, SP_R2 = SP2S2)

snv_v1 <- snv_v1 %>%
    select(-id) %>% 
    mutate_all(funs(norm = . / (sum(.) / 1000000))) %>% 
    bind_cols(select(snv_v1, id), .)

# proportion of cells in each bin (relative to pre-sort)
# DP1S2, DP2S2, INT1S2, INT2S2, PS_R1, PS_R2, SP1S2, SP2S2
bin_prop <- c(0.092, 0.091, 0.127, 0.122, 1, 1, 0.596, 0.603)

# multiply each bin count by bin proportion
snv_v1 <- bind_cols(select(snv_v1, header = id, DP_R1:SP_R2), 
                    data.frame(mapply(`*`, 
                                      select(snv_v1, DP_R1_norm:SP_R2_norm), 
                                      bin_prop, SIMPLIFY = FALSE)))


# SNV second sequencing run (v2), four bins
snv_v2 <- read.csv('../../processed_data/snv/snv_v2_all_alignments.csv')
snv_v2 <- snv_v2 %>% 
    select(Hi.R1:Lo.R2) %>% 
    mutate_all(funs(norm = . / (sum(., na.rm = T) / 1000000))) %>% 
    bind_cols(select(snv_v2, id), .)

# Hi.R1, Hi.R2, IntHi.R1, IntHi.R2, IntLo.R1, IntLo.R2, Lo.R1, Lo.R2
bin_prop <- c(0.070, 0.075, 0.029, 0.032, 0.032, 0.033, 0.227, 0.252)

# multiply each bin count by bin proportion
snv_v2 <- bind_cols(select(snv_v2, header = id, Hi.R1:Lo.R2), 
                    data.frame(mapply(`*`, 
                                      select(snv_v2, Hi.R1_norm:Lo.R2_norm), 
                                      bin_prop, SIMPLIFY = FALSE)))

###############################################################################
# Filtering
###############################################################################

hi_read_threshold <- 10
rep_agreement <- 0.20

# read filter
snv_v1 <- snv_v1 %>%
    mutate(v1_R1_sum = DP_R1 + INT_R1 + SP_R1,
           v1_R2_sum = DP_R2 + INT_R2 + SP_R2) %>%
    filter(v1_R1_sum >= hi_read_threshold, v1_R2_sum >= hi_read_threshold)

snv_v2 <- snv_v2 %>%
    mutate(v2_R1_sum = Hi.R1 + IntHi.R1 + IntLo.R1 + Lo.R1,
           v2_R2_sum = Hi.R2 + IntHi.R2 + IntLo.R2 + Lo.R2) %>% 
    filter(v2_R1_sum >= hi_read_threshold, v2_R2_sum >= hi_read_threshold) 

print(paste("Number of sequences after read filter (v1, v2):", 
            nrow(snv_v1), nrow(snv_v2)))

# index agreement before replicate filter
snv_v1 <- snv_v1 %>% 
    # calculate index
    mutate(v1_index_R1 = (DP_R1_norm * 0 + INT_R1_norm * 0.85 + 
                              SP_R1_norm * 1) / (DP_R1_norm + INT_R1_norm + SP_R1_norm),
           v1_index_R2 = (DP_R2_norm * 0 + INT_R2_norm * 0.85 + 
                              SP_R2_norm * 1) / (DP_R2_norm + INT_R2_norm + SP_R2_norm),
           v1_R1_norm = DP_R1_norm + INT_R1_norm + SP_R1_norm,
           v1_R2_norm = DP_R2_norm + INT_R2_norm + SP_R2_norm,
           v1_norm = v1_R1_norm + v1_R2_norm) %>%
    # rep agreement
    filter(abs(v1_index_R1 - v1_index_R2) <= rep_agreement)

# index agreement
snv_v2 <- snv_v2 %>% 
    mutate(v2_index_R1 = (Hi.R1_norm * 0 + IntHi.R1_norm * 0.80 + 
                              IntLo.R1_norm * 0.95 + Lo.R1_norm * 1) / 
               (Hi.R1_norm + IntHi.R1_norm + IntLo.R1_norm + Lo.R1_norm), 
           v2_index_R2 = (Hi.R2_norm * 0 + IntHi.R2_norm * 0.80 + 
                              IntLo.R2_norm * 0.95 + Lo.R2_norm * 1) / 
               (Hi.R2_norm + IntHi.R2_norm + IntLo.R2_norm + Lo.R2_norm),
           v2_R1_norm = Hi.R1_norm + 
               IntHi.R1_norm + IntLo.R1_norm + Lo.R1_norm,
           v2_R2_norm = Hi.R2_norm + 
               IntHi.R2_norm + IntLo.R2_norm + Lo.R2_norm,
           v2_norm = v2_R1_norm + v2_R2_norm) 

# correlation for v2 before index filter
snv_v2 <- snv_v2 %>% 
    mutate(v2_index_R1_binary = ifelse(v2_index_R1 >= 0.5, 1, 0),
           v2_index_R2_binary = ifelse(v2_index_R2 >= 0.5, 1, 0))

# tetrachoric correlation, 0.94
weightedCorr(snv_v2$v2_index_R1_binary, 
             snv_v2$v2_index_R2_binary, 
             method = "polychoric", 
             weights = snv_v2$v2_norm)

gg <- snv_v2 %>%
    mutate(rep_quality = ifelse(abs(v2_index_R1 - v2_index_R2) <= 0.20, 
                                'high', 'low')) %>% 
    ggplot(aes(v2_index_R1, v2_index_R2)) + 
    geom_point(alpha = 0.10, aes(color = rep_quality)) +
    scale_color_manual(values = c('black', 'darkgrey')) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    geom_abline(slope = 1, intercept = -0.20, linetype = 'dashed') +
    geom_abline(slope = 1, intercept = 0.20, linetype = 'dashed') +
    labs(x = 'inclusion index (v2 biological replicate 1)', 
         y = 'inclusion index (v2 biological replicate 2)') +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 16, vjust = -2), 
          axis.title.y = element_text(size = 16, vjust = +4),
          axis.text.x = element_text(size = 14, color = 'grey20'),
          axis.text.y = element_text(size = 14, color = 'grey20'),
          axis.ticks.x = element_line(color = 'grey50'),
          axis.ticks.y = element_line(color = 'grey50'),
          axis.line.x = element_line(color = 'grey50'),
          axis.line.y = element_line(color = 'grey50'),
          plot.margin = unit(c(2,2,3,3),"mm"))
    # annotate('text', x = 0.10, y = 0.90, parse = T,
    #          label = paste('italic(r) ==', signif(corr[1], 2)), size = 5)

ggsave(paste0('../../figs/supplement/SF6A_snv_v2_replicates', plot_format), 
       gg, width = 6, height = 6)

# rep agreement
snv_v2 <- snv_v2 %>% 
    filter(abs(v2_index_R1 - v2_index_R2) <= rep_agreement)

###############################################################################
# Join data (v1 and v2)
###############################################################################
data_all <- full_join(snv_v1, snv_v2, by = 'header') %>% 
    # small substitutions so separate will work easier
    mutate(header = gsub('strand= ', 'strand=', header),
           header = gsub('>', '', header),
           all_norm = v1_norm + v2_norm) %>%
    separate(header, into = c('id', 'chr', 'strand', 'length'), sep = ' ') %>% 
    select(-chr, -strand, -length)

data_all$v1_index <- rowMeans(select(data_all, v1_index_R1, v1_index_R2))
data_all$v2_index <- rowMeans(select(data_all, v2_index_R1, v2_index_R2))

# tetrachoric correlation, 0.99
data_all <- data_all %>%
    filter(!is.na(v1_index), !is.na(v2_index)) %>%
    mutate(v1_index_binary = ifelse(v1_index >= 0.5, 1, 0),
           v2_index_binary = ifelse(v2_index >= 0.5, 1, 0))

weightedCorr(data_all$v1_index_binary, 
             data_all$v2_index_binary, 
             method = "polychoric", 
             weights = data_all$all_norm)

# correlation between v1 and v2
corr <- wtd.cor(data_all$v1_index, data_all$v2_index, data_all$all_norm)
gg <- ggplot(data_all, aes(v1_index, v2_index)) + geom_point(alpha = 0.10) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(x = 'inclusion index (SNV library v1)', 
         y = 'inclusion index (SNV library v2)') +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 16, vjust = -2), 
          axis.title.y = element_text(size = 16, vjust = +4),
          axis.text.x = element_text(size = 14, color = 'grey20'),
          axis.text.y = element_text(size = 14, color = 'grey20'),
          axis.ticks.x = element_line(color = 'grey50'),
          axis.ticks.y = element_line(color = 'grey50'),
          axis.line.x = element_line(color = 'grey50'),
          axis.line.y = element_line(color = 'grey50'),
          plot.margin = unit(c(2,2,3,3),"mm"))+
    theme(legend.position = 'none')
    # annotate('text', x = 0.95, y = 0.10, parse = T,
    #          label = paste0('italic(r)==', signif(corr[1], 2)), size = 5) +
    # annotate('text', x = 0.96, y = 0.05, parse = T,
    #          label = paste('italic(p) < 10^-16'), size = 5)

ggsave(paste0('../../figs/supplement/SF6B_snv_v1_v2_replicates', plot_format), 
       gg, width = 6, height = 6)

# read in updated ref
ref <- read.table(paste0('../../ref/snv/',
                         'snv_ref_formatted_converted_original_seq.txt'), 
                  sep = '\t', header = T)

# combine with data
data_all <- left_join(data_all, ref, by = 'id') %>% 
    arrange(ensembl_id, sub_id) %>% 
    # update the category to either control, natural, or mutant
    mutate(category = ifelse(endsWith(id, '000'), 
                             'natural', 'mutant'), 
           category = ifelse(endsWith(id, 'BRK'), 
                             'control', category),
           category = ifelse(endsWith(id, 'SKP'), 
                             'skipped_exon', category),
           category = ifelse(startsWith(ensembl_id, 'RANDOM-EXON'), 
                             'random_exon', category)) %>%
    # get rid of misc. data columns
    select(id, ensembl_id:category, v1_R1_sum:v1_index_R2, v1_index, 
           v2_R1_sum:v2_index_R2, v2_index)

# replace NaN with NA
data_all[data_all == 'NaN'] <- NA

###############################################################################
# Filtering on natural exons
###############################################################################
data <- data_all %>% 
    group_by(ensembl_id) %>% 
    # mutant must have corresponding natural that passed previous filters
    filter(any(sub_id == '000')) %>% 
    # calculate dPSI
    mutate(nat_v1_index_R1 = v1_index_R1[sub_id == '000'],
           nat_v1_index_R2 = v1_index_R2[sub_id == '000'],
           nat_v2_index_R1 = v2_index_R1[sub_id == '000'], 
           nat_v2_index_R2 = v2_index_R2[sub_id == '000'],
           v1_dpsi_R1 = v1_index_R1 - nat_v1_index_R1,
           v1_dpsi_R2 = v1_index_R2 - nat_v1_index_R2,
           v2_dpsi_R1 = v2_index_R1 - nat_v2_index_R1,
           v2_dpsi_R2 = v2_index_R2 - nat_v2_index_R2,
           nat_v1_index = (nat_v1_index_R1 + nat_v1_index_R2) / 2, 
           nat_v2_index = (nat_v2_index_R1 + nat_v2_index_R2) / 2,
           nat_seq = original_seq[sub_id == '000']) %>% 
    filter(abs(nat_v2_index_R1 - nat_v2_index_R2) <= rep_agreement) %>%
    ungroup() 

# get averages between replicates
data <- data %>%
    rowwise() %>%
    mutate(v1_dpsi = mean(c(v1_dpsi_R1, v1_dpsi_R2)),
           v2_dpsi = mean(c(v2_dpsi_R1, v2_dpsi_R2)),
           delta_dpsi = abs(v1_dpsi - v2_dpsi) ) %>%
    ungroup()

dpsi_threshold <- -0.50

data <- data %>% 
    mutate(strong_lof = ifelse(v2_dpsi <= dpsi_threshold, T, F),
           strong_lof_v1 = ifelse(v1_dpsi <= -0.50 & nat_v1_index >= 0.50, T, F),
           strong_lof = ifelse(is.na(v2_dpsi), NA, strong_lof),
           strong_lof_v1 = ifelse(is.na(v1_dpsi), NA, strong_lof_v1))


# natural index > 0.50
data <- data %>% 
    filter(nat_v2_index >= 0.50)

data_final <- data %>% 
    filter(nat_v2_index >= 0.5, category == 'mutant')

# keep control categories: SKP and RANDOM-EXON 
data_other <- data_all %>% 
    group_by(ensembl_id) %>% 
    filter(any(sub_id == 'SKP') | (any(ensembl_id == 'RANDOM-EXON')) ) %>%
    ungroup() %>% 
    mutate(category = 
               case_when(.$sub_id == 'SKP' ~ 'skipped control',
                         .$ensembl_id == 'RANDOM-EXON' ~ 'random nucleotides'))
data_other <- bind_rows(data_other, filter(data, category == 'control'))

# plot controls
gg <- data_other %>% 
    bind_rows(filter(data, category == 'natural')) %>% 
    ggplot(aes(v2_index)) + geom_density(aes(fill = category), alpha = 0.5) +
    scale_fill_discrete(labels = c('broken SD/SA control', 
                                   'wild-type sequences',
                                   'random nucleotides', 
                                   'skipped control')) +
    labs(x = 'inclusion index') +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))  +
    theme(legend.position = c(0.40, 0.75),
          legend.title = element_blank(),
          axis.title.x = element_text(size = 20, vjust = -2), 
          axis.title.y = element_text(size = 20, vjust = +4),
          axis.text.x = element_text(size = 14, color = 'grey20'),
          axis.text.y = element_text(size = 14, color = 'grey20'),
          axis.ticks.x = element_line(color = 'grey50'),
          axis.ticks.y = element_line(color = 'grey50'),
          axis.line.x = element_line(color = 'grey50'),
          axis.line.y = element_line(color = 'grey50'),
          plot.margin = unit(c(2,2,3,3),"mm"))

ggsave(paste0('../../figs/supplement/SF7_snv_controls', plot_format), 
       gg, width = 4.5, height = 4)

write.table(data, '../../processed_data/snv/snv_data_clean.txt', 
            sep = '\t', row.names = F, quote = F)