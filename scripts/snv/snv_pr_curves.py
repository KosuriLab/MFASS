### Generate data for PR and ROC curves for various external models ###

import pandas as pd
import numpy as np


def precision_recall(df, var1, var2):
    # assume var1 is ground truth
    # true positive, both calls agree
    num_true_pos = float(len(df[(df[var1] == True) & (df[var2] == True)]))
    # false positives, called as positive by var2 but not var1
    num_false_pos = float(len(df[(df[var1] == False) & (df[var2] == True)]))
    # true negative, both var1 and var2 are negative
    num_true_neg = float(len(df[(df[var1] == False) & (df[var2] == False)]))
    # false negative, called as negative by var2 but positive by var1
    num_false_neg = float(len(df[(df[var1] == True) & (df[var2] == False)]))

    # precision, how many selected items are relevant
    if num_true_pos + num_false_pos == 0:
        precision = float('NaN')
    else:
        precision = (num_true_pos / (num_true_pos + num_false_pos)) * 100
    # recall/sensitivity, how many relevant items are selected
    if num_true_pos + num_false_neg == 0:
        recall = float('NaN')
    else:
        recall = (num_true_pos / (num_true_pos + num_false_neg)) * 100
    # specificity, ability to correctly detect negatives 
    specificity = (num_true_neg / (num_true_neg + num_false_pos)) * 100
    
    return [precision, recall, specificity]


def pr_score_curve(df, direction, score_width, truth_var, score_var, score_min=None, score_max=None):
    # include max as threshold
    if score_min is None:
        # score_min = np.nanmin(df[score_var]).iloc[0] # returns Series
        score_min = np.nanmin(df[score_var])
    if score_max is None:
        # score_max = np.nanmax(df[score_var]).iloc[0]
        score_max = np.nanmax(df[score_var])
    score_thresholds = np.arange(score_min, score_max + score_width, score_width)
    precisions = []
    recalls = []
    for threshold in score_thresholds:
        # same direction means lower scores are associated with less splicing defects
        if direction == 'same':
            df['score_var_strong_lof'] = np.where(df[score_var] >= threshold, True, False)
        # opposite means higher scores are associated with less splicing defects
        elif direction == 'opposite':
            df['score_var_strong_lof'] = np.where(df[score_var] <= threshold, True, False)
        else:
            raise Exception("please indicate direction of score assocation, same or opposite")
        df.loc[np.isnan(df[score_var]), 'score_var_strong_lof'] = float('NaN')
        precision, recall, specificity = precision_recall(df, truth_var, 'score_var_strong_lof')
        precisions.append(precision)
        recalls.append(recall)
    result = pd.DataFrame({'threshold' : score_thresholds, 
        'precision' : precisions, 
        'recall' : recalls})
    return result



def pr_threshold_curve(df, truth_var, prediction_var):
    # for methods that predict dPSI, change threshold where variant is called strong loss of function and compare
    # performance
    thresholds = np.arange(-1, 0, 0.001)
    precisions = []
    recalls = []
    for threshold in thresholds:
        df['prediction_var_strong_lof'] = np.where(df[prediction_var] <= threshold, True, False)
        df.loc[np.isnan(df[prediction_var]), 'prediction_var_strong_lof'] = float('NaN')
        precision, recall, specificity = precision_recall(df, truth_var, 'prediction_var_strong_lof')
        precisions.append(precision)
        recalls.append(recall)
    result = pd.DataFrame({'threshold' : thresholds, 
        'precision' : precisions, 
        'recall' : recalls})
    return result



def run_pr_methods(df, intron=False):
    # if intron variants only, don't run HAL (only exonic variants)
    fathmm_noncoding_pr_curve = pr_score_curve(df, 
        score_min=0, score_max=1, 
        score_width=0.01, truth_var='strong_lof', 
        score_var='noncoding_score', direction='same')
    fathmm_noncoding_pr_curve['method'] = ['fathmm_noncoding'] * len(fathmm_noncoding_pr_curve)

    fathmm_coding_pr_curve = pr_score_curve(df, 
        score_min=0, score_max=1, 
        score_width=0.01, truth_var='strong_lof', 
        score_var='coding_score', direction='same')
    fathmm_coding_pr_curve['method'] = ['fathmm_coding'] * len(fathmm_coding_pr_curve)
    
    cadd_pr_curve = pr_score_curve(df, 
        score_width = 0.5, 
        truth_var='strong_lof', score_var='cadd_score', 
        direction='same')
    cadd_pr_curve['method'] = ['cadd'] * len(cadd_pr_curve)
    
    fitcons_pr_curve = pr_score_curve(df, 
        score_min=0, score_max=1, 
        score_width = 0.01, truth_var='strong_lof', 
        score_var='fitCons_score', direction='opposite')
    fitcons_pr_curve['method'] = ['fitcons'] * len(fitcons_pr_curve)
    
    dann_pr_curve = pr_score_curve(df, 
        score_min=0, score_max=1, 
        score_width = 0.001, truth_var='strong_lof', 
        score_var='dann_score', direction='same')
    dann_pr_curve['method'] = ['dann'] * len(dann_pr_curve)
    
    linsight_pr_curve = pr_score_curve(df, 
        score_min=0, score_max=1, 
        score_width = 0.001, truth_var='strong_lof', 
        score_var='linsight_score', direction='opposite')
    linsight_pr_curve['method'] = ['linsight'] * len(linsight_pr_curve)

    phastcons_pr_curve = pr_score_curve(df, 
        score_min=0, score_max=1, 
        score_width = 0.001, truth_var='strong_lof', 
        score_var='mean_cons_score', direction='opposite')
    phastcons_pr_curve['method'] = ['phastCons'] * len(phastcons_pr_curve)

    phylop_pr_curve = pr_score_curve(df, 
        score_min=-10.829, score_max=9.873, 
        score_width = 0.01, truth_var='strong_lof', 
        score_var='phylop_score', direction='opposite')
    phylop_pr_curve['method'] = ['phyloP'] * len(phylop_pr_curve)
    
    spanr_pr_curve = pr_threshold_curve(df, 'strong_lof', 'spanr_dpsi')
    spanr_pr_curve['method'] = ['spanr'] * len(spanr_pr_curve)
    
    if not intron:
        hal_pr_curve = pr_threshold_curve(df, 'strong_lof', 'hal_dpsi')
        hal_pr_curve['method'] = ['hal'] * len(hal_pr_curve)
        pr_curve_info = pd.concat([fathmm_noncoding_pr_curve, 
            fathmm_coding_pr_curve, cadd_pr_curve, fitcons_pr_curve,
            dann_pr_curve, linsight_pr_curve, spanr_pr_curve, hal_pr_curve,
            phylop_pr_curve, phastcons_pr_curve, spanr_pr_curve])
    else:
        pr_curve_info = pd.concat([fathmm_noncoding_pr_curve, 
            fathmm_coding_pr_curve, cadd_pr_curve, fitcons_pr_curve,
            dann_pr_curve, linsight_pr_curve, phastcons_pr_curve,
            phylop_pr_curve, spanr_pr_curve])
    return pr_curve_info



if __name__ == '__main__':

    dpsi_threshold = -0.50

    snv_spanr = pd.read_table('../../processed_data/snv/snv_SPANR_scores_capped.txt', 
        sep='\t', header=0)
    snv_spanr.rename(index=str, inplace=True, 
        columns={'dpsi_spanr_capped' : 'spanr_dpsi', 
        'dpsi_max_tissue' : 'spanr_dpsi_uncapped'})
    snv_spanr['spanr_strong_lof'] = np.where(snv_spanr['spanr_dpsi'] <= dpsi_threshold, 
        True, False)


    snv_exon_vars = pd.read_table('../../processed_data/snv/snv_HAL_scores.txt', 
        sep='\t', header=0)
    snv_exon_vars.rename(index=str, inplace=True, columns={'DPSI_pred' : 'hal_dpsi'})
    snv_exon_vars['hal_strong_lof'] = np.where(snv_exon_vars['hal_dpsi'] <= dpsi_threshold, 
        True, False)

    data_annot = pd.read_table('../../processed_data/snv/snv_func_annot.txt', 
        sep='\t', header=0)
    # only keep those with dPSI values
    data_annot = data_annot[~np.isnan(data_annot.v2_dpsi)]


    data_all = pd.merge(data_annot[['id', 'v2_dpsi', 'category', 'strong_lof', 
        'label', 'mean_cons_score', 'consequence', 'cadd_score', 'noncoding_score', 
        'coding_score', 'fitCons_score', 'dann_score', 'linsight_score', 'phylop_score']],
                       snv_exon_vars[['id', 'hal_dpsi', 'hal_strong_lof']],
                       how='left', on='id')

    data_all = pd.merge(data_all, snv_spanr[['id', 'spanr_dpsi', 'spanr_strong_lof']], 
                       how='left', on='id')

    # only keep mutants 
    data_all = data_all[data_all.category == 'mutant']
    data_exon = data_all[data_all.label == 'exon']
    data_intron = data_all[data_all.label != 'exon']
    pr_curve_all = run_pr_methods(data_all)
    pr_curve_all.to_csv('../../processed_data/snv/snv_models_pr_curves_all.txt', 
        sep='\t', index=False)

    pr_curve_exon = run_pr_methods(data_exon)
    pr_curve_exon.to_csv('../../processed_data/snv/snv_models_pr_curves_exon.txt', 
        sep='\t', index=False)

    pr_curve_intron = run_pr_methods(data_intron, intron=True)
    pr_curve_intron.to_csv('../../processed_data/snv/snv_models_pr_curves_intron.txt', 
        sep='\t', index=False)


