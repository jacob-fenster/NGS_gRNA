"""
This is to find the genes associated with gRNAs that have significant positive and negative lfc
"""
import os, sys, ipdb
import pandas as pd
import numpy as np

data_dir = '/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis'
os.makedirs(f"{data_dir}/split_genes/", exist_ok=True)

experiments = ['LB_t2-vs-LB_t1.sgrna_summary.txt', 'gluc_t2-vs-gluc_t1.sgrna_summary.txt', 
               'ace_t2-vs-ace_t1.sgrna_summary.txt', '10pCA_t2-vs-10pCA_t1.sgrna_summary.txt',
               '50pCA_t2-vs-50pCA_t1.sgrna_summary.txt']

split_gene_list = pd.DataFrame(columns=['condition', 'gene'])
for experiment in experiments:
    data = pd.read_csv(f"{data_dir}/{experiment}", sep='\t', index_col=0)
    summary = pd.DataFrame(columns=data.columns)
    genes = set(data['Gene'].tolist())
    for gene in genes:
        # filter out counts that start less than 50
        gene_df = data[(data["Gene"] == gene) & (data["control_mean"] > 50)]
        # filter signifcant lfc and signifcant p value for pos and neg LFC
        pos_df = gene_df[(gene_df["LFC"] > 0.8) & (gene_df["p.high"] < 0.05)]
        neg_df = gene_df[(gene_df["LFC"] < -0.8) & (gene_df["p.low"] < 0.05)]
        if len(pos_df) > 0 and len(neg_df) > 0:
            summary = pd.concat([summary, pos_df, neg_df], ignore_index=False)
            if np.isnan(split_gene_list.index.max()):
                split_gene_list.loc[0] = {'condition': f"{experiment.split('.')[0]}", 'gene': gene}
            else:
                split_gene_list.loc[split_gene_list.index.max() + 1] = {'condition': f"{experiment.split('.')[0]}", 'gene': gene}
    summary.to_csv(f"{data_dir}/split_genes/{experiment.split('.')[0]}.tsv", sep='\t')
split_gene_list.to_csv(f"{data_dir}/split_genes/split_genes_vs_condition.tsv", sep='\t')

