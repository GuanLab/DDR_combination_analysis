import pandas as pd
import numpy as np
import itertools
from scipy.stats import kruskal, wilcoxon, mannwhitneyu
from scikit_posthocs import posthoc_dunn
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr, spearmanr
import plotly.express as px
import os
from tqdm import tqdm

df_scrambled_results = pd.read_csv('./tmp/biomarker_response_cor_scrambled_1000.tsv', sep = '\t')
df_results = pd.read_csv('./tmp/biomarker_response_cor.tsv', sep = '\t')

import math
all_fdr_p = {'biomarker':[], 'fdr':[], 'combination':[], 'score':[], 'biomarker_type':[]}
all_fdr = []
#cor_cutoff =[0.1, 0.2, 0.3, 0.35, 0.4]
#p_cutoff = [0.05, 0.01, 0.001, 0.0001, 0.00001]
#
#n_top_cutoff = [0.1,0.09,0.08,0.07,0.06, 0.05, 0.04,0.03,0.02,0.01, 0.005, 0.001]
# for c in cor_cutoff:
for comb in ['gamma-peposertib', 'M4076-berzosertib']:
    for biotype in ['exp', 'cnv', 'snv', 'lof', 'ddr', 'coh_pat', 'lof_pat']:
        for s in ['aoc', 'bliss']:
            #tmp = df_results[df_results['combination'] == comb]
            #tmp = tmp[tmp['biomarker_type'] == biotype]
            #tmp = tmp[tmp['score'] == s]
            tmp_scrambled = df_scrambled_results[df_scrambled_results['combination'] == comb]
            tmp_scrambled = tmp_scrambled[tmp_scrambled['biomarker_type'] == biotype]
            tmp_scrambled = tmp_scrambled[tmp_scrambled['score'] == s]
            for b in tmp_scrambled['biomarker'].unique():
                cut = 0
                all = tmp_scrambled[tmp_scrambled['biomarker']==b].sort_values(by='cor(score,biomarker)', ascending=False)['cor(score,biomarker)']
                # we define fdr as the number of elements in all that are larger/smaller than cut

                # number of elements in all that are larger than cut
                if all.mean() > 0:
                    fp = [i for i in all if i < 0]
                    fdr = len(fp)/len(all)
                else:
                    fp = [i for i in all if i > 0]
                    fdr = len(fp)/len(all)
                all_fdr_p['fdr'].append(fdr)
                all_fdr_p['biomarker'].append(b)
                all_fdr_p['combination'].append(comb)
                all_fdr_p['score'].append(s)
                all_fdr_p['biomarker_type'].append(biotype)

                # we define noise: biomarkers with less than 95% of correlation with positive score
                if (not math.isnan(cut))&(fdr < 0.01):
                    print(comb, s, b, fdr)

all_fdr_p = pd.DataFrame.from_dict(all_fdr_p)

df_results = df_results.merge(all_fdr_p, on=['biomarker', 'combination', 'score', 'biomarker_type'])
df_results.to_csv('./tmp/biomarker_response_cor_fdr.tsv', sep = '\t', index=False)

