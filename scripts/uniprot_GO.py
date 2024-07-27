import os, sys, ipdb
import pandas as pd 

def search_gene_names(df, term):
    # Search in the specified columns
    result = df[
             df['Gene Names (synonym)'].str.contains(term, na=False) | 
             df['Gene Names (ordered locus)'].str.contains(term, na=False) | 
             df['Gene Names (primary)'].str.contains(term, na=False)]
    if len(result) > 0:
        index = result.index[0]
        go = result.loc[index, 'Gene Ontology (GO)']
        go_process = result.loc[index, 'Gene Ontology (biological process)']
        go_component = result.loc[index, 'Gene Ontology (cellular component)']
        go_function = result.loc[index, 'Gene Ontology (molecular function)']
        return go, go_process, go_component, go_function
    else:
        return None, None, None, None
data_go = pd.read_csv("/Users/jfenster/Documents/NGS_gRNA/data/putida_proteome/uniprotkb_proteome_UP000000556_2024_07_25.tsv", 
                      sep='\t', index_col=0)
split_folder = "/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/split_genes"
split_df = pd.read_csv(f"{split_folder}/split_genes_vs_condition_kegg.tsv", sep='\t', index_col=0)
no_result = []
for index in split_df.index:
    gene = split_df.loc[index, 'gene']
    go, go_process, go_component, go_function = search_gene_names(data_go, gene)
    if go is None:
        no_result.append(gene)
    split_df.loc[index, 'Gene Ontology (GO)'] = go
    split_df.loc[index, 'Gene Ontology (biological process)'] = go_process
    split_df.loc[index, 'Gene Ontology (cellular component)'] = go_component
    split_df.loc[index, 'Gene Ontology (molecular function)'] = go_function
split_df.to_csv(f"{split_folder}/split_genes_vs_condition_kegg_go.tsv", sep='\t')
print(no_result)
breakpoint()
