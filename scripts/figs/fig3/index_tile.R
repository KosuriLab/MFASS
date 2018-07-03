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

fig_folder <- '../../../figs/fig2/'

###############################################################################
#### index, average index at each binned relative position, heatmap tile ###
### Figure 2B ###
###############################################################################

data <- data %>%
    mutate(rel_pos_binned = cut(rel_position_scaled, 
                                breaks = seq(-0.80, 1.80, 0.01)))

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

# # save legend separately
# g_legend <- function(a.gplot){
#     tmp <- ggplot_gtable(ggplot_build(a.gplot))
#     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#     legend <- tmp$grobs[[leg]]
#     legend
# }
# legend <- g_legend(index_tile_with_legend)
# tiff(paste0(fig_folder, 'legend_2C.tiff'), width = 25, height = 27.5, units = 'mm',
#      res = hi_res)
# grid.newpage()
# grid.draw(legend)
# dev.off()

index_tile <- index_tile_with_legend + theme(legend.position = 'none')

ggsave(paste0(fig_folder, 'snv_index_tile', plot_format), index_tile,
       width = 11, height = 1, units = 'in')
