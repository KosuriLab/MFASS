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

# GO term enrichment in SDV genes
ensembl <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
gene_info <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 
                                           'go_id', 'name_1006', 'definition_1006'), 
                            mart = ensembl,
                            filters = 'ensembl_gene_id',
                            values = data$ensembl_gene_id)
source("http://bioconductor.org/biocLite.R")
# biocLite('topGO')
# biocLite("org.Hs.eg.db")
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