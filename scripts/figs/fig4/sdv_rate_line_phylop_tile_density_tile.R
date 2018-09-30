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

###############################################################################
# Figure 4C
###############################################################################
### index, binned by relative position ###
data$label_renamed <- factor(data$label, 
                             levels=c("upstr_intron", "exon", "downstr_intron"), 
                             labels=c("5' intron", "exon", "3' intron"))

data <- data %>%
    mutate(rel_pos_binned = cut(rel_position_scaled, 
                                breaks = seq(-0.80, 1.80, 0.01)))

# color_intron = "#b90c0d"
# color_exon = "black"
# 
# group.colors <- c("5' intron" = color_intron,  
#                   "exon" = color_exon,  
#                   "3' intron" = color_intron) 

# SDVs by position
# absolute position relative to the center of the acceptor splice site at
# the end of the intron, 3' end intron and 5' end of exon
# acceptor splice site is at end of upstream intron and is 23bp long, extending
# 3bp into the exon. Center is -9 from intron/exon boundary
# data <- data %>% 
#     mutate(rel_pos_acc = ifelse(strand == '+',
#                                 rel_position - (intron1_len - 8),
#                                 rel_position - (intron2_len - 8)))

### SDV rate line graphs
tmp <- data %>% 
    filter(!is.na(rel_position_scaled), category != 'control') %>% 
    mutate(rel_pos_binned = cut(rel_position_scaled, breaks = seq(-0.80, 1.80, 0.02))) %>% 
    group_by(rel_pos_binned) %>% 
    summarise(sdv_rate = (length(which(strong_lof == T)) / n()) * 100)

tmp %>% 
    mutate(rel_pos_binned = as.numeric(rel_pos_binned)) %>% 
    ggplot(aes(rel_pos_binned, sdv_rate)) + 
    geom_line(aes(group = 1)) +
    # geom_area(aes(y = sdv_rate), fill = 'darkgrey') +
    scale_y_log10() + annotation_logticks(sides = 'l') +
    scale_x_discrete(breaks = levels(tmp$rel_pos_binned)[c(41, 91)],
                     labels = c(0, 1)) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = '', 
         y = '',
         color = '')

ggsave(paste0(fig_folder, 'sdv_rate_line_graph', plot_format),
       width = 11, height = 3, units = 'in')


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

ggsave(filename = paste0(fig_folder, "snv_phylop_tile", plot_format), 
       plot = phylop_tile,
       width = 11, height = 1, units = 'in')

# phylop_legend <- get_legend(phylop_tile)
# ggsave(paste0(fig_folder, 'phylop_legend', plot_format),
#        phylop_legend, width = 1, height = 3, units = 'in')

###############################################################################
### SNV density tile ###
###############################################################################
ref <- read.table('../../../ref/snv/snv_ref_formatted_converted.txt',
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
          # legend.position = 'none',
          axis.line = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0,0,0,0),'in')) +
    viridis::scale_fill_viridis() +
    labs(x = '', y = '', fill = '')

ggsave(paste0(fig_folder, 'snv_density_tile', plot_format), snp_density_tile,
       width = 11, height = 1, units = 'in')
