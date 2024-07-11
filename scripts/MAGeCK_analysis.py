import os, sys, glob, pdb, time
import pandas as pd

file = '/Users/jacobfenster/Documents/NGS_gRNA/output/20231126_exact_counts/MAGeCK_test/LB_t0_t2_Rep1_mle-noday0_perm2.gene_summary.txt'
df = pd.read_csv(file, sep='\t')
mle_columns = ['Gene', 'sgRNA', 'LB_t1_vs_t0_R1|beta', 'LB_t1_vs_t0_R1|z', 'LB_t1_vs_t0_R1|p-value', 'LB_t1_vs_t0_R1|fdr', 'LB_t1_vs_t0_R1|wald-p-value', 'LB_t1_vs_t0_R1|wald-fdr']
