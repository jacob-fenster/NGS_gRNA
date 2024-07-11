import sys, pdb
import pandas as pd
#this one extracts the seeds from each 
   

gRNA_file_path = "/Users/jacobfenster/Documents/NGS_gRNA/data/PPutida_gRNA_Libraries/Ordered_filtered_databases/Lib3_dCas9_filtered_ordered.csv"
df = pd.read_csv(gRNA_file_path, index_col=0)
names, oligos = df.index.tolist(), df['gRNA'].tolist()
seeds = [oligo[65:85] for oligo in oligos]
genes = [name.split(':')[0] for name in names]
seed_df = pd.DataFrame(data={'gene':genes, 'seed':seeds}, index=names)
seed_df.index.name = 'sgRNA'
seed_df.to_csv("/Users/jacobfenster/Documents/NGS_gRNA/databases/Lib3_dCas9_78K_seeds.csv")