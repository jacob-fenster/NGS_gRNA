from bioservices.kegg import KEGG
import pandas as pd
import os, glob

k = KEGG()
k.organism = 'ppu'

split_folder = "/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/split_genes"
experiments = ['LB_t2-vs-LB_t1.tsv', 'gluc_t2-vs-gluc_t1.tsv', 
               'ace_t2-vs-ace_t1.tsv', '10pCA_t2-vs-10pCA_t1.tsv',
               '50pCA_t2-vs-50pCA_t1.tsv']
split_data_file = "/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/split_genes/split_genes_vs_condition.tsv"
split_data = pd.read_csv(split_data_file, sep='\t', index_col=0)
for index in split_data.index:
    gene = split_data.loc[index, 'gene']
    data = k.get(f"{k.organism}:{gene}")
    if type(data) != str:
        print(f"error with {gene}: {data}")
    else:
        gene_dict = k.parse(data)
        try:
            split_data.loc[index, 'NAME'] = gene_dict['NAME']
        except:
            continue
        try:
            split_data.loc[index, 'ORTHOLOGY'] = str(gene_dict['ORTHOLOGY'])
        except:
            continue
        try:
            pathway_id = k.get_pathway_by_gene(gene, k.organism)
            split_data.loc[index, 'PATHWAY ID'] = str(pathway_id)
        except:
            continue
split_data.to_csv(f"{split_folder}/split_genes_vs_condition_kegg.tsv", sep='\t')

        
        