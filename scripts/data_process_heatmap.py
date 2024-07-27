import os, sys, glob, ipdb
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

data_dir = '/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis'
os.makedirs(f"{data_dir}/clustering", exist_ok=True)
experiments = ['LB_t2-vs-LB_t1.gene_summary.txt', 'gluc_t2-vs-gluc_t1.gene_summary.txt', 
               'ace_t2-vs-ace_t1.gene_summary.txt', '10pCA_t2-vs-10pCA_t1.gene_summary.txt',
               '50pCA_t2-vs-50pCA_t1.gene_summary.txt']
exp_tags = ['LB', 'glucose', 'acetate', 'low_pCA', 'high_pCA']
pval_tr = 0.05
test_data = pd.read_csv(f"{data_dir}/{experiments[0]}", sep='\t', index_col=0)
master_df = pd.DataFrame(index=test_data.index)
for experiment, tag in zip(experiments, exp_tags):
    data = pd.read_csv(f"{data_dir}/{experiment}", sep='\t', index_col=0)
    pos_data = data[data["pos|lfc"] > 0]
    neg_data = data[data["neg|lfc"] <= 0]
    
    # set the lfc values equal to zero where p values are below tr
    pos_data.loc[pos_data['pos|p-value'] > pval_tr, 'pos|lfc'] = 0
    neg_data.loc[neg_data['neg|p-value'] > pval_tr, 'neg|lfc'] = 0
    # combine the datasets 
    pos_data = pos_data['pos|lfc']
    neg_data = neg_data['neg|lfc']
    processed_data = pd.concat([pos_data, neg_data])
    master_df[f"lfc|{tag}"] = processed_data
master_df.to_csv(f"{data_dir}/clustering/20240722_gene_lfc.tsv", sep='\t')

# now normalize 
breakpoint()
scaler = MinMaxScaler(feature_range=(-1, 1)) 
master_df.iloc[:, :] = scaler.fit_transform(master_df.iloc[:, :])
master_df.to_csv(f"{data_dir}/clustering/20240722_gene_lfc_minmaxscale.tsv", sep='\t')
breakpoint()