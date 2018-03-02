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
hi_res <- 600

data <- read.table('../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
  filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)

###############################################################################
# Figure 3A, ExAC global allele frequency 
###############################################################################
data <- data %>%
  mutate(AF_bin = cut(AF, 
                      breaks = c(0, 0.00001, 0.000025, 0.0005, 0.005, 1), 
                      include.lowest = T,
                      labels = c('Singleton', 'AC = 2-3', 
                                 'AC = 4-10', '0.05-0.5%', '>0.5%' ))) 

sdv_mutant_v2 <- data %>% 
  filter(nat_v2_index >= 0.50, category == 'mutant') %>%
  mutate(sdv = ifelse(v2_dpsi <= -0.50, T, F)) %>% 
  group_by(category, AF_bin) %>% 
  summarise(num_sdv = length(which(sdv == T)),
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
