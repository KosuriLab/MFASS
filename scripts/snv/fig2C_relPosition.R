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

hi_res <- 300
plot_format_main <- '.png'
plot_format <- '.png'

data <- read.table('../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

dpsi_threshold <- -0.50

###############################################################################
# Figure 2C
###############################################################################
### index, binned by relative position for boxplot ###
data$label_renamed <- factor(data$label, 
                             levels=c("upstr_intron", "exon", "downstr_intron"), 
                             labels=c("5' intron", "exon", "3' intron"))

data <- data %>%
    mutate(rel_pos_binned = cut(rel_position_scaled, 
                                breaks = seq(-0.80, 1.80, 0.01)))

color_intron = "#b90c0d"
color_exon = "black"

group.colors <- c("5' intron" = color_intron,  
                  "exon" = color_exon,  
                  "3' intron" = color_intron) 

# SDVs by position
# absolute position relative to the center of the acceptor splice site at
# the end of the intron, 3' end intron and 5' end of exon
# acceptor splice site is at end of upstream intron and is 23bp long, extending
# 3bp into the exon. Center is -9 from intron/exon boundary
data <- data %>% 
  mutate(rel_pos_acc = ifelse(strand == '+',
                              rel_position - (intron1_len - 8),
                              rel_position - (intron2_len - 8)))

### SDV rate line graphs
tmp <- data %>% 
  filter(!is.na(rel_position_scaled), category != 'control') %>% 
  mutate(rel_pos_binned = cut(rel_position_scaled, breaks = seq(-0.80, 1.80, 0.02))) %>% 
  group_by(rel_pos_binned) %>% 
  summarise(sdv_rate = (length(which(strong_lof == T)) / n()) * 100)

tmp %>% 
  ggplot(aes(rel_pos_binned, sdv_rate)) + 
  geom_point() +
  geom_line(aes(group = 1)) +
  geom_area() +
  scale_y_log10() + annotation_logticks(sides = 'l') +
  scale_x_discrete(breaks = levels(tmp$rel_pos_binned)[c(41, 91)],
                   labels = c(0, 1)) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(colour = "grey20", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
        # axis.line.x = element_blank(),
        # axis.ticks.x = element_blank()) +
  labs(x = '', 
       y = '',
       color = '')

ggsave(paste0('../../figs/snv/fig2C_sdv_rate_line_graph', '.png'),
       width = 11, height = 3, units = 'in')

###############################################################################
#### index, average index at each binned relative position, heatmap tile ###
###############################################################################
index_tile_with_legend <- data %>%
    filter(category == "mutant") %>%
    filter(!is.na(rel_pos_binned)) %>% 
    group_by(rel_pos_binned) %>% 
    summarise(mean_dpsi_per_rel_pos = mean(v2_dpsi, na.rm = T)) %>% 
    ggplot(aes(x = rel_pos_binned, y = 0.5)) + 
    geom_tile(aes(fill = mean_dpsi_per_rel_pos)) +
    viridis::scale_fill_viridis(option = "viridis", direction = -1, 
                                breaks = seq(-1, 0.1, 0.25), 
                                limits = c(-1, 0.10)) +
    labs(x = '', y = '', fill = '') +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.line = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "in"))

# save legend separately
g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}
legend <- g_legend(index_tile_with_legend)
tiff('../../figs/snv/legend_2C.tiff', width = 25, height = 27.5, units = 'mm', 
     res = hi_res)
grid.newpage()
grid.draw(legend)
dev.off()

index_tile <- index_tile_with_legend + theme(legend.position = 'none')

ggsave(paste0('../../figs/snv/fig2C_snv_index_tile', '.png'),
    width = 11, height = 1, units = 'in')


###############################################################################
### phastCons scores for all natural sequences in library
###############################################################################
create_positions <- function(ensembl_id, chr, start, end, strand) {
    all_positions = data.frame()
    if(strand == '-'){ counter = 170}
    else { counter = 1}
    for(i in seq(start, end + 1)) {
        if(counter > 170 | counter <= 0) { next }
        id = paste0(ensembl_id, '_pos', counter)
        if(strand == '-') {counter = counter - 1}
        else { counter = counter + 1}
        all_positions <- bind_rows(all_positions,
                                   data.frame(chr = chr, start = i, 
                                              end = i + 1, id = id))
    }
    return(all_positions)
}

# for each natural exon background, generate data frame with each genomic position
# of the intron-exon-intron construct
if(!file.exists('../../processed_data/snv/snv_nat_cons_scaled.rds')) {
    data %>% 
        filter(category == 'natural') %>% 
        group_by(ensembl_id) %>% 
        do(data.frame(create_positions(.$ensembl_id, .$chr, .$start_hg38_0based, 
                                       .$end_hg38_0based, .$strand))) %>% 
        ungroup() %>% 
        select(-ensembl_id) %>% 
        write.table('../../processed_data/snv/snv_nat_positions.bed',
                    sep = '\t', quote = F, row.names = F)
    
    # grab conservation at each base
    system(paste('bash',
                 '../run_phastCons.sh',
                 '../../processed_data/snv/snv_nat_positions.bed', 
                 '../../processed_data/snv/snv_nat_cons_scores_all.bed'))
    
    # read in data
    snv_nat_cons <- read.table('../../processed_data/snv/snv_nat_cons_scores_all.bed', 
                                sep = '\t', header = F,
                                col.names = c('name', 'size', 'bases_covered', 
                                              'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>% 
        separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = F) %>%
        mutate(rel_position = as.numeric(gsub('pos', '', position))) %>% 
        select(-(position:mean0), id = name) %>% 
        # join other information necessary to calculate relative position
        left_join(select(data, ensembl_id, intron1_len, 
                         exon_len, intron2_len, strand) %>% 
                      distinct(), 
                  by = 'ensembl_id')
    
    rel_position <- function(intron1_len, exon_len, intron2_len, strand, rel_position) {
        in_interval <- function(start, end, x) {
            if (start <= x & x <= end) { return(TRUE) }
            else {return(FALSE)}
        }
        if (strand == '-' | strand == -1) {
            # set appropriate lengths
            upstr_intron_len <- intron2_len
            downstr_intron_len <- intron1_len
        }
        else{
            upstr_intron_len <- intron1_len
            downstr_intron_len <- intron2_len
        }
        
        regions <- data.frame(label = c('upstr_intron', 'exon', 'downstr_intron'),
                              start = c(1, upstr_intron_len + 1, upstr_intron_len + exon_len + 1),
                              end = c(upstr_intron_len, upstr_intron_len + exon_len, 
                                      upstr_intron_len + exon_len + downstr_intron_len))
        # select region the SNP falls in
        label <- regions %>% 
            rowwise() %>% 
            filter(in_interval(start, end, x = rel_position))
        
        start <- label$start[1]
        end <- label$end[1]
        label <- label$label[1]
        
        # get distance from end of upstream intron/exon boundary
        boundary <- upstr_intron_len
        distance <- rel_position - boundary
        
        # normalize to feature length
        if (label == 'downstr_intron') {
            scaled_distance <- 1 + (distance - exon_len) / (end - start)
        }
        else {
            scaled_distance <- distance / (end - start + 1)
        }
        
        # additionally, find distance from respective intron/exon boundary
        # left side of boundary is negative, right side is positive
        if (label == 'upstr_intron') {
            rel_pos_feature <- rel_position - end - 1
        }
        if (label == 'exon') {
            # closer to downstream intron
            if ( rel_position - start + 1 >= end - rel_position + 1) {
                # negative position, left side of boundary
                rel_pos_feature <- rel_position - end - 1 
            }
            else {
                # closer to upstream intron, right side of boundary, positive
                rel_pos_feature <- rel_position - start + 1
            }
        }
        if (label == 'downstr_intron') {
            rel_pos_feature <- rel_position - start + 1
        }
        return(paste(label, rel_pos_feature, scaled_distance, sep = ':'))
    }
    
    snv_nat_cons <- snv_nat_cons %>% 
        select(intron1_len, exon_len, intron2_len, strand, rel_position, id) %>%
        na.omit() %>%
        rowwise() %>%
        mutate(rel_position_info = rel_position(intron1_len, exon_len, 
                                                intron2_len, strand, rel_position)) %>%
        separate(rel_position_info, 
                 c('label', 'rel_position_feature', 'rel_position_scaled'), 
                 sep = ':', convert = T) %>%
        select(id, label, rel_position_feature, rel_position_scaled) %>%
        left_join(snv_nat_cons, ., by = 'id') %>% distinct()
    
    snv_nat_cons <- snv_nat_cons %>% 
        mutate(rel_pos_binned = cut(rel_position_scaled, 
                                    breaks = seq(-.80, 1.80, 0.01)))
    # save as RDS so factor is retained
    saveRDS(snv_nat_cons, '../../processed_data/snv/snv_nat_cons_scaled.rds')
} else {
    snv_nat_cons <- readRDS('../../processed_data/snv/snv_nat_cons_scaled.rds')
}

###############################################################################
# phyloP tile
###############################################################################
phylop_tile <- data %>% 
    filter(!is.na(rel_pos_binned), category == 'mutant') %>% 
    group_by(rel_pos_binned) %>% 
    summarise(mean_phylop = mean(phylop_score)) %>% 
    ggplot(aes(x = rel_pos_binned, y = 0.5)) + geom_tile(aes(fill = mean_phylop)) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = 'none',
          axis.line = element_blank(),
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "in")
    ) +
    labs(y = '', x = '', fill = '') +
    viridis::scale_fill_viridis()
phylop_tile
ggsave(filename = paste0("../../figs/snv/fig2C_snv_phylop_tile", '.png'), 
       plot = phylop_tile,
       width = 11, height = 1, units = 'in')

###############################################################################
### SNP density tile ###
###############################################################################
ref <- read.table('../../ref/snv/snv_ref_formatted_converted.txt',
                  sep = '\t', header = T)

total_snps <- ref %>%
    filter(sub_id != '000', sub_id != 'BRK') %>%
    nrow()

snp_density_tile <- ref %>%
    mutate(rel_pos_binned = cut(rel_position_scaled, breaks = seq(-.80, 1.80, 0.01))) %>%
    filter(!is.na(rel_pos_binned), sub_id != '000', sub_id != 'BRK') %>%
    group_by(rel_pos_binned) %>%
    summarise(snp_density = n() / total_snps) %>%
    ggplot(aes(rel_pos_binned, 0.5)) + geom_tile(aes(fill = snp_density)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position = 'none',
          axis.line = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0,0,0,0),'in')) +
    viridis::scale_fill_viridis() +
    labs(x = '', y = '', fill = '')

ggsave(paste0('../../figs/snv/fig2C_snv_density', '.png'),
    width = 11, height = 1, units = 'in')
