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

data <- read.table('../../processed_data/snv/snv_func_annot.txt',
                   sep = '\t', header = T)
data <- data %>% 
    filter(category == 'mutant', (!is.na(v2_dpsi)), nat_v2_index >= 0.5)


# Clinvar
sdvs <- filter(data, strong_lof == T)
# hg19 coordinates
clinvar <- read.table('../../ref/clinvar_variant_summary.txt',
                      header = F, sep = '\t', fill = T)
names <- c('AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol',	'HGNC_ID',
           'ClinicalSignificance',	'ClinSigSimple', 'LastEvaluated', 
           'rsid', 'nsv_esv_dbvar',	'RCVaccession',	'PhenotypeIDS',	'PhenotypeList',
           'Origin', 'OriginSimple', 'Assembly', 'ChromosomeAccession',
           'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele',
           'Cytogenetic', 'ReviewStatus', 'NumberSubmitters', 'Guidelines',
           'TestedInGTR', 'OtherIDs', 'SubmitterCategories', 'VariationID')
names(clinvar) <- names
# join based on chromosome and SNP position
clinvar_sdvs <- clinvar %>% 
    mutate(Start = as.numeric(Start)) %>% 
    semi_join(sdvs, by = c('Chromosome' = 'chr', 'Start' = 'snp_position')) %>% 
    left_join(select(sdvs, id, v2_dpsi, strong_lof, consequence, chr, snp_position, ref_allele, alt_allele),
              by = c('Chromosome' = 'chr', 'Start' = 'snp_position')) %>% 
    select(id, v2_dpsi, consequence, Type, GeneSymbol, ClinicalSignificance,
           Chromosome, Start, strong_lof,
           ref_allele, ReferenceAllele, alt_allele, AlternateAllele) %>% 
    filter(ref_allele == ReferenceAllele)

clinvar_snvs <- clinvar %>% 
    mutate(Start = as.numeric(Start)) %>% 
    semi_join(data, by = c('Chromosome' = 'chr', 'Start' = 'snp_position')) %>% 
    left_join(select(data, id, v2_dpsi, strong_lof, consequence, chr, snp_position, ref_allele, alt_allele),
              by = c('Chromosome' = 'chr', 'Start' = 'snp_position')) %>% 
    select(id, v2_dpsi, consequence, Type, GeneSymbol, ClinicalSignificance,
           Chromosome, Start, strong_lof,
           ref_allele, ReferenceAllele, alt_allele, AlternateAllele) %>% 
    filter(ref_allele == ReferenceAllele)


df1 <- clinvar_snvs %>% 
    filter(strong_lof == T) %>% 
    group_by(ClinicalSignificance) %>% 
    tally() %>% 
    mutate(pct = n / sum(n),
           sdv = 'SDV (n = 8)')

df2 <- clinvar_snvs %>% 
    filter(strong_lof == F) %>% 
    group_by(ClinicalSignificance) %>% 
    tally() %>% 
    mutate(pct = n / sum(n),
           sdv = 'all SNVs (n = 141)')

df3 <-  clinvar_snvs %>% 
    filter(strong_lof == F) %>% 
    group_by(ClinicalSignificance) %>% 
    tally() %>% 
    mutate(pct = n / sum(n),
           sdv = 'nonSDV (n = 133)')

bind_rows(df1, df2, df3) %>% 
    ggplot(aes(ClinicalSignificance, pct * 100)) + 
    geom_bar(stat = 'identity', position = 'dodge', aes(fill = sdv)) + 
    scale_fill_manual(values = c('grey40', 'black', 'darkblue')) +
    labs(y = 'percentage', x = 'Clinvar significance') +
    # theme(axis.text.x = element_text(angle = 45)) + 
    coord_flip()

ggsave('../../figs/clinvar_sdvs.png',
       width = 8, height = 4, units = 'in')
