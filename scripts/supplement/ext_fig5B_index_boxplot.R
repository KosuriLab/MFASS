load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'forcats', 'gridExtra', 
          'grid', 'Unicode')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

hi_res <- 300
plot_format_main <- '.tiff'
plot_format <- '.tiff'

data <- read.table('../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

dpsi_threshold <- -0.50

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

index_boxplot <- data %>% 
    filter(category == "mutant") %>% 
    filter(!is.na(rel_pos_binned)) %>% 
    ggplot(aes(rel_pos_binned, v2_dpsi)) + 
    geom_boxplot(aes(color = label_renamed), outlier.size = 0.10, 
                 outlier.colour = "black", notch = FALSE) +
    theme(legend.position = "none",
          axis.text.x = element_blank(), 
          axis.title.y = element_text(margin = 
                                        margin(0, 0, -65, -65), size = 20),
          axis.title.x =  element_text(margin = 
                                         margin(0, 0, -45, -65), size = 20),
          axis.text.y = element_text(colour = "grey20", size = 16),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = '', 
         y = expression(paste(Delta, ' inclusion index')),
         color = '') + 
    scale_color_manual(values = group.colors) + ylim(-1, 1)

ggsave(paste0('../../figs/snv/ext_fig5B_index_boxplot', plot_format),
       width = 11, height = 4, units = 'in')
