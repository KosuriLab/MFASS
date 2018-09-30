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

plot_format <- '.png'
hi_res <- 600

data <- read.table('../../processed_data/snv/snv_func_annot.txt', sep = '\t', header = T)
data <- data %>%
    filter(category == 'mutant', nat_v2_index >= 0.5, !(is.na(v2_dpsi)))

data <- data %>% 
    mutate(sdv = ifelse(nat_v2_index >= 0.5 & v2_dpsi <=-0.5, 'SDV', 'non-SDV'))

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

# get genes in terms
sig_terms <- result$GO.ID[1:4]
go_term_genes <- genesInTerm(GOdata, sig_terms)
go_df <- data.frame(GO.ID = rep.int(names(go_term_genes), lengths(go_term_genes)),
                    ensembl_gene_id = unlist(go_term_genes))
go_df <- go_df %>% 
    left_join(dplyr::select(result, GO.ID, Term))

go_df <- go_df %>% 
    left_join(dplyr::select(data, ensembl_gene_id, external_gene_name) %>% distinct())

gene_exon_info <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'ensembl_exon_id',
                                                'rank'), 
                                 mart = ensembl,
                                 filters = 'ensembl_gene_id',
                                 values = data$ensembl_gene_id)

gene_exon_num <- gene_exon_info %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(total_exons  = max(rank))

go_summary <- data %>% 
    filter(ensembl_gene_id %in% go_df$ensembl_gene_id) %>% 
    mutate(sdv = ifelse(v2_dpsi <= -0.50, 'SDV', 'non-SDV')) %>% 
    group_by(ensembl_gene_id, sdv) %>% 
    tally() %>%
    mutate(pct = (n / sum(n)) * 100) %>% 
    left_join(go_df, by = 'ensembl_gene_id') %>% 
    left_join(gene_exon_num, by = 'ensembl_gene_id')

go_summary %>% 
    distinct(GO.ID, ensembl_gene_id, total_exons) %>% 
    group_by(GO.ID) %>% 
    summarise(mean_exon_num = mean(total_exons, na.rm = T),
              median = median(total_exons, na.rm = T)) %>% 
    left_join(dplyr::select(go_df, GO.ID, Term))

data <- data %>% 
    left_join(gene_exon_num, by = 'ensembl_gene_id')

sdv_info <- data %>% 
    filter(category == 'mutant', strong_lof == T) %>% 
    dplyr::select(id, chr, snp_position, snp_position_hg38_1based,
                  ref_allele, alt_allele, strand, label,
                  ensembl_gene_id, external_gene_name, total_exons,
                  snp_id, exon_id_new,
                  v2_dpsi, nat_v2_index,
                  AF_gnomad, consequence, mean_phastCons_score,
                  phylop_score, SIFT)

write.table(sdv_info, '../../processed_data/snv/sdv_info.txt',
            sep = '\t', quote = F, row.names = F)
