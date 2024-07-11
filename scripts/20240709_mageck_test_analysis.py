"""
This is an initial test to look at the MAGeCK 'test' function data in a "tabular volcano plot format"
last run 20240709
"""
import os, sys, ipdb, glob
import pandas as pd
import numpy

output_dir = "output/20240709_tabular_volcano"
gene_summary_files = glob.glob(f"output/20231226_main_analysis/*.gene_summary.txt")
for gene_summary_file in gene_summary_files:
    df = pd.read_csv(gene_summary_file, index_col=0, sep='\t')
    filtered_df = df[(df['neg|p-value'] <= 0.05) | (df['pos|p-value'] <= 0.5)]
    # these are off the cuff cutoffs. Need to logically pick the values 
    filtered_df = filtered_df[(filtered_df['neg|lfc'] <= -1) | (filtered_df['pos|lfc'] >= 0.5)]
    pos_sig = filtered_df[(filtered_df['pos|p-value'] <= 0.5) & (filtered_df['pos|lfc'] >= 0.5)]
    neg_sig = filtered_df[(filtered_df['neg|p-value'] <= 0.5) & (filtered_df['neg|lfc'] <= -1)]
    pos_sig = pos_sig.sort_values(by="pos|lfc", ascending=False)
    neg_sig = neg_sig.sort_values(by="pos|lfc", ascending=True)
    pos_sig.to_csv(f"{output_dir}/positive/{os.path.basename(gene_summary_file).split('.txt')[0]}_pos_sig.csv")
    neg_sig.to_csv(f"{output_dir}/negative/{os.path.basename(gene_summary_file).split('.txt')[0]}_neg_sig.csv")

