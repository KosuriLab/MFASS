load_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
  if(length(new_pkgs)) install.packages(new_pkgs)
  for(pkg in pkgs){
    suppressWarnings(suppressMessages(library(pkg, character.only = T)))
  }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'grid', 'gtable', 'ggsignif')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.tiff'
hi_res <- 300

data <- read.table('../../processed_data/snv/snv_func_annot.txt', sep = '\t', header = T)
data <- data %>%
  filter(category == 'mutant', nat_v2_index >= 0.5, !(is.na(v2_dpsi)))

##############################
# Figure 3b: pLI
##############################

dpsi_threshold <- -0.50

num_intolerant_lof <- data %>% 
  filter(category == 'mutant', pLI_high == TRUE, v2_dpsi <= dpsi_threshold ) %>% nrow()
num_intolerant_not_lof <- data %>% 
  filter(category == 'mutant', pLI_high == TRUE, v2_dpsi > dpsi_threshold ) %>% nrow()
num_tolerant_lof <- data %>%
  filter(category == 'mutant', pLI_high == FALSE, v2_dpsi <= dpsi_threshold ) %>% nrow()
num_tolerant_not_lof <- data %>% 
  filter(category == 'mutant', pLI_high == FALSE, v2_dpsi > dpsi_threshold ) %>% nrow()

df <- matrix(c(num_intolerant_lof, num_tolerant_lof,
               num_intolerant_not_lof, num_tolerant_not_lof), 
             nrow = 2, 
             dimnames = list(c('pLI >= 0.90', 'pLI < 0.90'),
                             c('dPSI <= -0.50', 'dPSI > -0.50')))

intolerant_percent =  num_intolerant_lof / (num_intolerant_lof + num_intolerant_not_lof) * 100
tolerant_percent = num_tolerant_lof / (num_tolerant_not_lof + num_tolerant_lof) * 100
percent.df <-data.frame( fraction_of_strong_LoF_genes = c(intolerant_percent, tolerant_percent), 
                         tolerance = c('intolerant', 'tolerant'))

color_pli = "black"

percent.df %>%
  ggplot(aes(tolerance, fraction_of_strong_LoF_genes)) + 
  geom_col(width = 0.65, color = color_pli, fill = color_pli) + 
  ylab("% SDV") + xlab("pLI") + ylim(0,5) +
  geom_hline(yintercept = 1050/27733*100, linetype = "dashed", color = "grey40") +
  scale_y_continuous(expand = c(0, 0)) + 
  expand_limits(y = 4.8) +
  # coord_equal(1/1.125) +
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50"),
        axis.title.y = element_text(size = 16, color = "black", vjust = 2),
        axis.title.x = element_text(size = 20, color = "black", vjust = -0.5),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "grey20"),
        axis.text.x = element_text(size = 13, color = "black")) 

fisher.test(df, alternative = 'two.sided')

ggsave(paste0("../../figs/snv/fig3B_pLI_enrichment", plot_format), 
       width = 2.25, height = 3.5, units = 'in', dpi = hi_res)
