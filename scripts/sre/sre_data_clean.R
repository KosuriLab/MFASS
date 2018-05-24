### Process raw data with appropriate normalization and filters ###

load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

ref <- read.table('../../ref/sre/sre_ref_formatted_converted.txt', sep = '\t',
                  header = T, colClasses = c('sub_id' = 'character'))

smn1 <- read.csv('../../processed_data/sre/smn1/smn1_all_alignments.csv') %>% 
    rename(DP_R1 = R1.DP, INT_R1 = R1.INT, PS_R1 = R1.PRESORT, SP_R1 = R1.SP, 
           DP_R2 = R2.DP, INT_R2 = R2.INT, PS_R2 = R2.PRESORT, SP_R2 = R2.SP)

dhfr <- read.csv('../../processed_data/sre/dhfr/dhfr_all_alignments.csv') %>% 
    rename(DP_R1 = R1.DP, INT_R1 = R1.INT, PS_R1 = R1.PRESORT, SP_R1 = R1.SP, DN_R1 = R1.DN,
           DP_R2 = R2.DP, INT_R2 = R2.INT, PS_R2 = R2.PRESORT, SP_R2 = R2.SP, DN_R2 = R2.DN)


# normalize to number of million reads in each sample
smn1 <- smn1 %>%
    select(-id) %>% 
    mutate_all(funs(norm = . / (sum(.) / 1000000))) %>% 
    bind_cols(select(smn1, id), .)

dhfr <- dhfr %>% 
    select(-id) %>% 
    mutate_all(funs(norm = . / (sum(.) / 1000000))) %>% 
    bind_cols(select(dhfr, id), .)

# proportion of cells sorted into each bin
smn1_bin_prop <- c(0.320, 0.057, 1, 0.443, 0.318, 0.057, 1, 0.430)
dhfr_bin_prop <- c(0.024, 0.346, 0.077, 1, 0.363, 0.027, 0.365, 0.078, 1, 0.388)

# multiply each bin count by bin proportion
smn1 <- bind_cols(select(smn1, header = id, DP_R1:SP_R2), 
                  data.frame(mapply(`*`, select(smn1, DP_R1_norm:SP_R2_norm), 
                                    smn1_bin_prop, SIMPLIFY = FALSE)))

dhfr <- bind_cols(select(dhfr, header = id, DN_R1:SP_R2),
                  data.frame(mapply(`*`, select(dhfr, DN_R1_norm:SP_R2_norm),
                                    dhfr_bin_prop, SIMPLIFY = FALSE)))

# read sum across bins
smn1 <- smn1 %>% 
    mutate(R1_sum = DP_R1 + INT_R1 + SP_R1,
           R2_sum = DP_R2 + INT_R2 + SP_R2)

dhfr <- dhfr %>% 
    mutate(R1_sum = DP_R1 + INT_R1 + SP_R1,
           R2_sum = DP_R2 + INT_R2 + SP_R2)

# filter for low read count
low_read_threshold <- 5
smn1 <- smn1 %>% 
    filter(R1_sum >= low_read_threshold, R2_sum >= low_read_threshold)

dhfr <- dhfr %>% 
    filter(R1_sum >= low_read_threshold, R2_sum >= low_read_threshold)

# combine
data <- full_join(smn1, dhfr, by = 'header', suffix = c('_smn1', '_dhfr'))

# join reference 
data <- data %>% 
    rowwise %>% 
    mutate(id = unlist(strsplit(header, split = ' '))[1],
           id = ifelse(grepl('_', id), id, paste(id, '000', sep = '_'))) %>% 
    select(-header) %>% 
    left_join(ref, by = 'id') %>% 
    # reorder
    select(id:end_hg38, DP_R1_smn1:SP_R2_norm_dhfr) %>% 
    ungroup() %>% 
    arrange(ensembl_id, sub_id)

# convert blanks to NA
data[data == ''] <- NA

data <- data %>%
    mutate(index_R1_smn1 = (DP_R1_norm_smn1*0 + INT_R1_norm_smn1*0.85 + SP_R1_norm_smn1*1) / 
               (DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1),
           index_R2_smn1 = (DP_R2_norm_smn1*0 + INT_R2_norm_smn1*0.85 + SP_R2_norm_smn1*1) / 
               (DP_R2_norm_smn1 + INT_R2_norm_smn1 + SP_R2_norm_smn1),
           index_R1_dhfr = (DP_R1_norm_dhfr*0 + INT_R1_norm_dhfr*0.85 + SP_R1_norm_dhfr*1) / 
               (DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr),
           index_R2_dhfr = (DP_R2_norm_dhfr*0 + INT_R2_norm_dhfr*0.85 + SP_R2_norm_dhfr*1) / 
               (DP_R2_norm_dhfr + INT_R2_norm_dhfr + SP_R2_norm_dhfr),
           R1_norm_sum_smn1 = DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1,
           R2_norm_sum_smn1 = DP_R2_norm_smn1 + INT_R2_norm_smn1 + SP_R2_norm_smn1,
           norm_sum_smn1 = R1_norm_sum_smn1 + R2_norm_sum_smn1,
           R1_norm_sum_dhfr = DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr,
           R2_norm_sum_dhfr = DP_R2_norm_dhfr + INT_R2_norm_dhfr + SP_R2_norm_dhfr,
           norm_sum_dhfr = R1_norm_sum_dhfr + R2_norm_sum_dhfr)

data$index_smn1 <- rowMeans(select(data, index_R1_smn1, index_R2_smn1))
data$index_dhfr <- rowMeans(select(data, index_R1_dhfr, index_R2_dhfr))

# replace NaN with NA
data[data == 'NaN'] <- NA

rep_agreement <- 0.30
data <- data %>% 
    mutate(rep_quality = ifelse(abs(index_smn1 - index_dhfr) <= rep_agreement, 'high', 'low'),
           replicability_dhfr = ifelse(abs(index_R1_dhfr - index_R2_dhfr) <= rep_agreement, 'high', 'low'),
           replicability_smn1 = ifelse(abs(index_R1_smn1 - index_R2_smn1) <= rep_agreement, 'high', 'low'))

# calculate difference in splicing index between mutant and natural
data <- data %>% 
    group_by(ensembl_id) %>% 
    # filter out anything that doesn't have a corresponding natural sequence
    filter(any(sub_id == '000')) %>% 
    mutate(nat_index_R1_smn1 = index_R1_smn1[sub_id == '000'],
           nat_index_R2_smn1 = index_R2_smn1[sub_id == '000'],
           nat_index_R1_dhfr = index_R1_dhfr[sub_id == '000'], 
           nat_index_R2_dhfr = index_R2_dhfr[sub_id == '000'],
           dpsi_R1_smn1 = index_R1_smn1 - nat_index_R1_smn1,
           dpsi_R2_smn1 = index_R2_smn1 - nat_index_R2_smn1,
           dpsi_R1_dhfr = index_R1_dhfr - nat_index_R1_dhfr,
           dpsi_R2_dhfr = index_R2_dhfr - nat_index_R2_dhfr) %>%
    ungroup()

data <- data %>% 
    rowwise() %>% 
    mutate(dpsi_smn1 = mean(c(dpsi_R1_smn1, dpsi_R2_smn1)),
           dpsi_dhfr = mean(c(dpsi_R1_dhfr, dpsi_R2_dhfr)),
           nat_index_smn1 = mean(c(nat_index_R1_smn1, nat_index_R2_smn1)),
           nat_index_dhfr = mean(c(nat_index_R1_dhfr, nat_index_R2_dhfr))) %>% 
    ungroup()

write.table(data, '../../processed_data/sre/sre_data_clean.txt', sep ='\t', row.names = F)

data %>% 
    ggplot(aes(index_R1_smn1, index_R2_smn1)) + 
    geom_point(alpha = 0.10) +
    geom_density2d(color = 'black') + 
    geom_abline(linetype = 'dashed', intercept = 0.30) +
    geom_abline(linetype = 'dashed', intercept = -0.30) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
    scale_color_manual(values = c('black', '#0033CC')) +
    labs(x = 'inclusion index (Rep. 1)', y = 'inclusion index (Rep. 2)') +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 20, vjust = -2), 
          axis.title.y = element_text(size = 20, vjust = +4),
          axis.text.x = element_text(size = 14, color = 'grey20'),
          axis.text.y = element_text(size = 14, color = 'grey20'),
          axis.ticks.x = element_line(color = 'grey50'),
          axis.ticks.y = element_line(color = 'grey50'),
          axis.line.x = element_line(color = 'grey50'),
          axis.line.y = element_line(color = 'grey50'),
          plot.margin = unit(c(2,2,3,3),"mm"))

# overall SDV rate
data %>% 
    filter(rep_quality == 'high', seq_type == 'mut') %>% 
    mutate(sdv = ifelse((dpsi_smn1 <= -0.50) & nat_index_smn1 >= 0.5, T, F)) %>% 
    group_by(sdv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n))

data %>% 
    filter(rep_quality == 'high', seq_type == 'mut') %>% 
    mutate(sdv = ifelse((dpsi_dhfr <= -0.50) & nat_index_dhfr >= 0.5, T, F)) %>% 
    group_by(sdv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n))
