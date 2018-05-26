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

data <- read.table('../../processed_data/snv/snv_func_annot.txt', sep = '\t', header = T)
data <- data %>%
  filter(category == 'mutant', nat_v2_index >= 0.5, !(is.na(v2_dpsi)))


##############################
# Figure 3b: pLI
##############################

# num_intolerant_lof <- data %>% 
#   filter(category == 'mutant', pLI_high == TRUE, v2_dpsi <= dpsi_threshold ) %>% nrow()
# num_intolerant_not_lof <- data %>% 
#   filter(category == 'mutant', pLI_high == TRUE, v2_dpsi > dpsi_threshold ) %>% nrow()
# num_tolerant_lof <- data %>%
#   filter(category == 'mutant', pLI_high == FALSE, v2_dpsi <= dpsi_threshold ) %>% nrow()
# num_tolerant_not_lof <- data %>% 
#   filter(category == 'mutant', pLI_high == FALSE, v2_dpsi > dpsi_threshold ) %>% nrow()

dpsi_threshold <- -0.50
data <- data %>% 
    mutate(sdv = ifelse(v2_dpsi <= dpsi_threshold, 'SDV', 'non-SDV'))

counts <- data %>% 
    filter(!is.na(pLI_high)) %>% 
    mutate(tolerance = ifelse(pLI_high == T, 'intolerant', 'tolerant')) %>% 
    group_by(tolerance, sdv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n))

color_pli = "black"
counts %>%
    filter(sdv == 'SDV') %>%
    ggplot(aes(tolerance, pct*100)) +
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


# current matrix arrangement tests for significantly greater proportion of non-SDVs
# in intolerant genes, equivalent to significantly fewer proportion of SDVS in intolerant
# genes if you rearranged the matrix
df <- matrix(counts$n, nrow = 2)
fisher.test(df, alternative = 'two.sided')

ggsave(paste0("../../figs/snv/fig3B_pLI_enrichment", plot_format), 
       width = 2.25, height = 3.5, units = 'in', dpi = hi_res)


# compare rate of PTVs in pLI genes to SDVs
data <- data %>% 
    mutate(ptv = ifelse(consequence %in% c('stop_gained', 
                                           'splice_acceptor_variant',
                                           'splice_donor_variant'),
                        'PTV', 'non-PTV'))

counts <- data %>% 
    filter(!is.na(pLI_high)) %>% 
    mutate(tolerance = ifelse(pLI_high == T, 'intolerant', 'tolerant')) %>% 
    group_by(tolerance, ptv) %>% 
    tally() %>% 
    mutate(pct = n / sum(n)) %>% 
    rename(type = ptv) %>% 
    bind_rows(rename(counts, type = sdv))
    
counts %>% 
    filter(type == 'PTV' | type == 'SDV') %>% 
    ggplot(aes(tolerance, pct*100)) + 
    geom_bar(stat = 'identity', position = 'dodge', aes(fill = type)) + 
    scale_fill_manual(values = c('black', 'darkblue')) + 
    labs(x = 'pLI', y = 'percentage', fill = 'variant type') + ylim(0,5) +
    geom_hline(yintercept = 1050/27733*100, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 340/27733*100, linetype = "dashed", color = "blue") +
    scale_y_continuous(expand = c(0, 0)) + 
    expand_limits(y = 4.8) +
    # coord_equal(1/1.125) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(size = 16, color = "black", vjust = 2),
          axis.title.x = element_text(size = 20, color = "black", vjust = -0.5),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "grey20"),
          axis.text.x = element_text(size = 13, color = "black")) 


# GO term enrichment in SDV genes
ensembl <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
gene_info <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                                           'go_id', 'name_1006', 'definition_1006'), 
                            mart = ensembl,
                            filters = 'ensembl_gene_id',
                            values = data$ensembl_gene_id)
source("http://bioconductor.org/biocLite.R")
biocLite('topGO')
biocLite("org.Hs.eg.db")
library(topGO)

nonsdv_genes <- data %>% filter(sdv == 'non-SDV') %>% .$ensembl_gene_id
sdv_genes <- data %>% filter(sdv == 'SDV') %>% .$ensembl_gene_id
all_genes <- data$ensembl_gene_id
gene_list <- factor(as.integer(all_genes %in% sdv_genes))
names(gene_list) <- all_genes

GOdata <- new('topGOdata',
              ontology = 'BP',
              allGenes = gene_list,
              annot = annFUN.org,
              mapping = 'org.Hs.eg.db',
              ID = 'ensembl',
              nodeSize = 10)
# over-representation of GO terms within the group of SDV genes, each category
# tested independently
go_fisher <- runTest(GOdata, algorithm = 'classic', statistic = 'fisher')
# significant genes correspond to genes of interest
# non-trivial nodes are GO categories which have at least one significant gene annotated
result <- GenTable(GOdata, classicFisher = go_fisher)
