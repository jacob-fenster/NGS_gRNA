import os
import sys
import pandas as pd

def count_fasta_entries_and_get_basenames(fasta_file_list):
    entry_counts = []
    basenames = []
    for fasta_file in fasta_file_list:
        # Initialize count for the current file
        count = 0
        # Open and read FASTA file to count entries
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    count += 1   
        entry_counts.append(count)
        # Get base name of the file without extension
        base_name = os.path.basename(fasta_file)
        base_name_no_ext = os.path.splitext(base_name)[0]
        basenames.append(base_name_no_ext)
    return entry_counts, basenames

output_file, fasta_file_list = sys.argv[1], sys.argv[2:]
entry_counts, basenames = count_fasta_entries_and_get_basenames(fasta_file_list)
df = pd.DataFrame({'Nseq_80': entry_counts}, index=basenames)
df.to_csv(output_file)