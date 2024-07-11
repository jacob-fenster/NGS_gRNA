import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
import pdb

#this script filters sequences that are below a %ID and above a coverage tr
def read_blast_tsv(file_path):
#inputs a .tsv file from blastp (-outfmt "7 delim=    ) and outputs a pandas df
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Find the line starting with '# Fields:' to get column names
    columns = []
    for line in lines:
        if line.startswith("# Fields:"):
            # Remove '# Fields:' prefix and then split by ', '
            columns = line[len("# Fields: "):].strip().split(', ')
            break
    # Remove lines starting with '#' for actual data
    data_lines = [line.strip() for line in lines if not line.startswith('#')]
    #convert to list format for dataframe assembly
    rows = [dict(zip(columns, data_line.split('\t'))) for data_line in data_lines]
    #convert to correct type
    df = pd.DataFrame(rows)
    return df

def filter_fasta_by_dataframe(filtered_df, fasta_path, output_basename):
    # Load the FASTA file
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    # Create an empty list to hold filtered sequences
    filtered_sequences = []
    # Set of unique "subject acc.ver" strings from the dataframe
    subject_acc_set = set(filtered_df['subject acc.ver'].apply(str))
    for fasta in fasta_sequences:
        description = fasta.description
        # Check if part of the sequence ID is in "subject acc.ver" set
        if any(substring in description for substring in subject_acc_set):
            filtered_sequences.append(fasta)
        else:
            pdb.set_trace()
            print(1)
    # Write the filtered sequences to a new FASTA file
    output_path = f"{output_basename}_filtered.fasta"
    SeqIO.write(filtered_sequences, output_path, "fasta")
    return

#def filter_fasta_from_BLASTp(filtered_df, fasta, output_file):
#this one accepts a filtered dataframe blastp output, the fasta that made that output, and the outputfile and created a filtered fasta file

    #return 

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 blastp-tsv_filter.py <input_fasta_file> <input_tsv_file> <%ID tr> <coverage tr>")
        sys.exit(1)

fasta_path, tsv_path, ID_filter, coverage_filter = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
output_basename = os.path.basename(fasta_path).split('.')[0]

#convert .tsv into dataframe
tsv = read_blast_tsv(tsv_path)
#keep the first row for addition after filtering
df2 = pd.DataFrame(data=None, columns=tsv.columns)
df2.loc[0] = tsv.loc[0]
#filter dataframe and add back the query identity row
tsv['% identity'] = tsv['% identity'].astype(float)
tsv['% query coverage per subject'] = tsv['% query coverage per subject'].astype(float)
filtered_tsv = tsv[(tsv['% identity']<=ID_filter) & (tsv['% query coverage per subject']>=coverage_filter)]
filtered_tsv = pd.concat([df2, filtered_tsv], ignore_index=True)
filtered_tsv.to_csv(output_basename+"_filtered.tsv", sep='\t')
#filter the fasta file given the filtered dataframe
filter_fasta_by_dataframe(filtered_tsv, fasta_path, output_basename)