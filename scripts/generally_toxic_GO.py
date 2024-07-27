import os, sys, ipdb
import pandas as pd 

import seaborn as sns
import matplotlib.pyplot as plt

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
toxic = pd.read_csv("/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/clustering/generally_toxic.tsv", 
                    sep='\t', index_col=0)
no_result = []
for index in toxic.index:
    gene = index 
    go, go_process, go_component, go_function = search_gene_names(data_go, gene)
    if go is None:
        no_result.append(gene)
    toxic.loc[index, 'Gene Ontology (GO)'] = go
    toxic.loc[index, 'Gene Ontology (biological process)'] = go_process
    toxic.loc[index, 'Gene Ontology (cellular component)'] = go_component
    toxic.loc[index, 'Gene Ontology (molecular function)'] = go_function
toxic.to_csv("/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/clustering/generally_toxic_GO.tsv", sep='\t')
# sort into categories 
go_df = pd.DataFrame(columns=['Gene Ontology (biological process)'])
for index in toxic.index:
    try:
        go_list = toxic.loc[index, 'Gene Ontology (biological process)'].split(';')
    except AttributeError:
        go_list = ['None']
    for entry in go_list:
        if len(go_df.index) == 0:
            go_df.loc[0, 'Gene Ontology (biological process)'] = entry
        else:
            go_df.loc[len(go_df.index)+1, 'Gene Ontology (biological process)'] = entry

go_counts = go_df['Gene Ontology (biological process)'].value_counts()
go_counts = go_counts.drop('None')
go_desc = go_counts.index.tolist()
new_desc = []
for desc in go_desc:
    new = desc.split(' [GO:')[0]
    new_desc.append(new)
go_counts.index = new_desc
breakpoint()
# Plot the top 20 GO terms
top_n = 10
plt.figure(figsize=(5, 5))
go_counts.head(top_n).plot(kind='bar')
plt.title(f'Top {top_n} GO Term Counts')
plt.xlabel('GO Terms')
plt.ylabel('Counts')
plt.xticks(rotation=90)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()