import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
import pdb 

def fasta_MSA_summary_clusatlo(MSA_file):
#this inputs the path to a clustal-W MSA file and reads with BioPython
#generating a df with description versus summary stats
#assumes that first entry is the query
    names, seq, ID, cov, gap_ratio = [], [], [], [], []
    with open(MSA_file, 'r') as handle: #read the input fasta to lists 
        for record in SeqIO.parse(handle, "fasta"):
            names.append(record.description)
            seq.append(str(record.seq))
    full_len = len(seq[0]) #length of each entry
    query_seq = seq[0] #assumes first entry is query sequence
    query_len = float(len(query_seq.replace('-',''))) #length of query no gaps
    for i in range(0, len(names)): #iterate over MSAs, i's
        identical, align_len, aligned, gap = 0, 0, 0, 0 #identical are identical algined columns, not gaps. align_length is length of alignment getting rid of common gaps. Aligned are identical+nonidentical columns. Gaps are AA in query but gaps in target, doesn't count insertions
        for c in range(0, full_len): #iterate over MSA columns, c's
            if query_seq[c].isalpha() and seq[i][c].isalpha(): #test to see if both columns are alphanumerical
                if query_seq[c] == seq[i][c]: #test if identical
                    identical += 1
                aligned += 1
                align_len += 1
            elif (query_seq[c].isalpha() and seq[i][c]=='-') or (query_seq[c]=='-' and seq[i][c].isalpha()):
                align_len += 1
            if (query_seq[c].isalpha() and seq[i][c]=='-'):
                gap += 1
        ID.append(100*(identical/align_len))
        cov.append(100*(aligned/query_len))
        gap_ratio.append(100*(gap/query_len))
    df = pd.DataFrame(data={'Descriptions':names, '%ID':ID, 'coverage':cov, 'gap_ratio':gap_ratio})
    return df 

def filter_MSA_df(df, ID_tr, cov_tr, gap_ratio_tr):
#this inputs a datafame of summary stats from a MSA and list of filter paramters
#filter_list is ['column]
    #hold query row
    #pdb.set_trace()
    df2 = pd.DataFrame(data=None, columns=df.columns)
    df2.loc[0] = df.loc[0]
    filtered_df = df[(df['%ID']<=ID_tr) & (df['coverage']>=cov_tr) & (df['gap_ratio']>=gap_ratio_tr)]
    filtered_df = pd.concat([df2, filtered_df], ignore_index=True) #add back the query column
    return filtered_df

def filter_fasta_by_df(filtered_df, fasta_path, output_basename):
#this filters the given fasta file by the 'Descriptions' column of the input dataframe
    # Load the FASTA file
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    # Create an empty list to hold filtered sequences
    filtered_sequences = []
    # Set of unique "subject acc.ver" strings from the dataframe
    subject_acc_set = set(filtered_df['Descriptions'].apply(str))
    for fasta in fasta_sequences:
        description = fasta.description
        # Check if part of the sequence ID is in "subject acc.ver" set
        if any(substring in description for substring in subject_acc_set):
            filtered_sequences.append(fasta)
    # Write the filtered sequences to a new FASTA file
    output_path = f"{output_basename}_filtered.fasta"
    SeqIO.write(filtered_sequences, output_path, "fasta")
    return

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 clustalo_fasta_summary.py <input_MSAfasta_file> <%ID tr> <coverage tr> <gap ratio tr>")
        sys.exit(1)
MSA_file, ID_tr, cov_tr, gap_ratio_tr = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
if '10239' in MSA_file: #check to see if orthoDB ASFV level. specific to my nonsense naming scheme
    output_basename = os.path.basename(MSA_file).split('__')[0].split('_')[1] + os.path.basename(MSA_file).split('.')[0].split('10239')[1].split('_')[0]
else: #catch all for other file types
    output_basename = os.path.basename(MSA_file).split('.')[0].split('_rearranged_clustalo')[0]
full_basename = os.path.basename(MSA_file).split('.')[0]
df = fasta_MSA_summary_clusatlo(MSA_file)
df.to_csv(output_basename+'.csv')
filtered_df = filter_MSA_df(df, ID_tr, cov_tr, gap_ratio_tr)
filtered_df.to_csv(output_basename+'_filtered.csv')
filter_fasta_by_df(filtered_df, MSA_file, full_basename)
