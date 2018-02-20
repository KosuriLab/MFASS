readRenviron("~/.Renviron")
setwd(Sys.getenv("SPLICEVAR"))

### Graph PR and ROC curves for various external models ###

load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'grid', 'pROC')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

hi_res <- 300

color_cadd = '#FF9912'
color_dann = '#941494'
color_fathmm = '#108C44'
color_fitcons = '#6AA5CD'
color_linsight = '#ABB9B9'
color_spanr =  '#ED1E24'
color_hal = '#000080'
color_phastcons <- 'black'
color_phylop <- 'turquoise'

## ROC curves ##

data <- read.table('../../processed_data/exac/exac_func_annot.txt',
                   sep = '\t', header = T)
spanr <- read.table('../../processed_data/exac/exac_SPANR_scores_capped.txt',
                    sep = '\t', header = T)
hal <- read.table('../../processed_data/exac/exac_HAL_scores.txt',
                  sep = '\t', header = T)

data <- data %>% 
    left_join(dplyr::select(spanr, id, dpsi_spanr_capped), by = 'id') %>% 
    left_join(dplyr::select(hal, id, DPSI_pred), by = 'id') %>% 
    filter(category != 'control', nat_v2_index >= 0.5, (!is.na(v2_dpsi))) %>% 
    mutate(sdv = ifelse(strong_lof == T, 1, 0)) %>% 
    filter(!is.na(sdv))

run_roc <- function(df, type) {
    roc_objs <- list()
    if(type == 'exon'){
        roc_objs[['hal']] <- roc(df$sdv, df$DPSI_pred, ci = T)
    }
    roc_objs[['spanr']] <- roc(df$sdv, df$dpsi_spanr_capped, ci = T)
    roc_objs[['dann']] <- roc(df$sdv, df$dann_score, ci = T)
    roc_objs[['linsight']] <- roc(df$sdv, df$linsight_score, ci = T)
    roc_objs[['fitcons']] <- roc(df$sdv, df$fitCons_score, ci = T)
    roc_objs[['fathmm']] <- roc(df$sdv, df$coding_score, ci = T)
    roc_objs[['cadd']] <- roc(df$sdv, df$cadd_score, ci = T)
    roc_objs[['phastCons']] <-roc(data$sdv, data$mean_cons_score, ci = T)
    roc_objs[['phylop']] <- roc(data$sdv, data$phylop_score, ci = T)
    return(roc_objs)
}

all_roc_objs <- run_roc(data, type = 'all')
intron_roc_objs <- run_roc(filter(data, label != 'exon'), type = 'intron')
exon_roc_objs <- run_roc(filter(data, label == 'exon'), type = 'exon')

plot_rocs <- function(roc_objs, exon = F) {
    plot(roc_objs[['spanr']], col = color_spanr, legacy.axes = T, xlab = '', ylab = '')  # changed
    plot(roc_objs[['dann']], add = T, col = color_dann)
    plot(roc_objs[['linsight']], add = T, col = color_linsight)
    plot(roc_objs[['fitcons']], add = T, col = color_fitcons)
    plot(roc_objs[['fathmm']], add = T, col = color_fathmm)
    plot(roc_objs[['cadd']], add = T, col = color_cadd)
    plot(roc_objs[['phastCons']], add = T, col = color_phastcons)
    plot(roc_objs[['phylop']], add = T, col = color_phylop)
    if(exon){
        plot(roc_objs[['hal']], add = T, col = color_hal)
    }
}

methods <- c('HAL', 'SPANR', 'DANN', 'LINSIGHT', 'fitCons', 'FATHMM-MKL', 'CADD')
colors <- c(color_hal, color_spanr, color_dann, color_linsight, color_fitcons,
            color_fathmm, color_cadd)

tiff(paste0("../../figs/exac/roc_all.tiff"), width = 3.34, height = 3.44, units = 'in', res = 300)
plot_rocs(all_roc_objs, exon = F)
# legend("bottomright", legend = methods[-1], col = colors[-1], lty = 1, cex = 0.5)
dev.off()

tiff(paste0("../../figs/exac/roc_intron.tiff"), width = 3.34, height = 3.44, units = 'in', res = 300)
plot_rocs(intron_roc_objs, exon = F)
# legend("bottomright", legend = methods, col = colors, lty = 1, cex = 0.5)
dev.off()

tiff(paste0("../../figs/exac/roc_exon.tiff"), width = 3.34, height = 3.44, units = 'in', res = 300)
plot_rocs(exon_roc_objs, exon = T)
# legend("bottomright", legend = methods, col = colors, lty = 1, cex = 0.5)
dev.off()

# ROC for phyloP and phastCons
roc_phastCons <- roc(data$sdv, data$mean_cons_score)
roc_phylop <- roc(data$sdv, data$phylop_score)
plot(roc_phastCons, col = 'red', legacy.axes = T)
plot(roc_phylop, add = T, col = 'black', legacy.axes = T)
legend('bottomright', legend = c('phastCons', 'phyloP'), col = c('red', 'black'),
       lty = 1, cex = 0.5)
