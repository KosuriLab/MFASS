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

color_unknown = '#FF6103' 
color_probably_damaging = '#002F6C'
color_possibly_damaging = '#6AA5CD'
color_benign = '#9D1309'
color_noannot = '#AAAAAA'

################################
# PolyPhen and SIFT annotations
################################

data <- data %>% 
  separate(SIFT, into = c('SIFT', 'SIFT_other'), sep = '[(]', remove = T) %>% 
  separate(PolyPhen, into = c('PolyPhen', 'PolyPhen_other'), sep = '[(]', remove = T) %>% 
  dplyr::select(-SIFT_other, -PolyPhen_other)

################
#### Polyphen
################

data %>%
  filter(consequence == "missense_variant", v2_dpsi <= -0.5, !is.na(PolyPhen)) %>%
  ggplot(aes(x = PolyPhen, fill = PolyPhen)) + 
  geom_histogram(stat="count", width = 0.6) + 
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 180) +
  labs(x = '', #PolyPhen prediction\n for missense SDVs
       y = 'count') +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50"),
        axis.title.y = element_text(vjust = 20, size = 20),
        axis.title.x = element_text(vjust = -20, size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
  ) +
  scale_fill_manual(values = c(color_unknown, 
                               color_probably_damaging, 
                               color_possibly_damaging, 
                               color_benign,
                               color_noannot)) +
  scale_x_discrete(labels = c('no annotation', 
                              'benign', 
                              'possibly\ndamaging', 
                              'probably\ndamaging', 
                              'unknown')) +
  coord_flip()

ggsave(paste0('../../figs/snv/fig4A_polyphen_SDVs_missense', plot_format), 
       width = 5, height = 3, dpi = hi_res) 



############
#### SIFT
############

color_tolerated =  '#FF6103'
color_deleterious = '#9D1309'
color_noannot = '#AAAAAA'

data %>%
  filter(consequence == "missense_variant", v2_dpsi <= -0.5, !is.na(SIFT)) %>%
  filter(SIFT != "tolerated_low_confidence", SIFT != "deleterious_low_confidence") %>%
  ggplot(aes(x = SIFT, fill = SIFT)) +
  geom_histogram(stat="count", width = 0.6) + 
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 180) +
  # geom_jitter(alpha = 0.1) +
  # stat_summary(fun.y = median, geom = "point", size = 3) +
  labs(x = '',  #SIFT prediction\nfor missense SDVs
       y = 'count') +
  scale_x_discrete(labels = c('no annotation', 
                              'deleterious', 
                              'tolerated')) +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 20, vjust = 60),
        axis.title.x = element_text(vjust = -20, size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14)
  ) +
  coord_flip() +
  scale_fill_manual(values = c(color_tolerated, 
                               color_deleterious, 
                               color_noannot))

ggsave(paste0('../../figs/snv/fig4B_SIFT_SDVs_missense', plot_format), 
       width = 5, height = 2.1, dpi = hi_res) 


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
    ggplot(aes(ClinicalSignificance, pct)) + 
    geom_bar(stat = 'identity', position = 'dodge', aes(fill = sdv)) + 
    labs(y = 'percentage', x = 'Clinvar significance') +
    # theme(axis.text.x = element_text(angle = 45)) + 
    coord_flip()
