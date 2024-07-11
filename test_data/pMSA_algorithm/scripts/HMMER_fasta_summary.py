import os
import sys
import numpy as np
import pandas as pd
import pdb 
# INCOMPLETE
#this script inputs a fasta formatted MSA from HMMER .sto output
#first entry must be query 
def fasta_MSA_summary(MSA_file):
    query_name, query_seq = "", ""
    query_found = True
    entries, sequence, header = {}, "", None
    # Read the FASTA file and populate the 'entries' dictionary. Use first entry as query
    with open(MSA_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries[header] = sequence
                    if query_found: #check if the first line has been found
                        query_name, query_seq = header, sequence
                        query_found = False
                header, sequence = line, ""
            else:
                sequence += line
        if header:
            entries[header] = sequence
    pdb.set_trace()
    # Calculate stats on entries
    coverage, ID = '', ''
    query_len = float(len(query_seq.replace('.', '').replace('-','')))
    for header, sequence in entries.items():
        coverage_seq = ''.join(c for c in sequence if not c.islower()).replace('.', '').replace('-','')
        coverage_len = float(len(coverage_seq))
        coverage = coverage_len / query_len #this calculation might overestimate because it looks at alignment rows, not matches
        alignment_len = sequence.replace('.', '').replace('-','')
        matches = 0
        for i in range(0, len(sequence)):
            if sequence[i] == query_seq[i] and sequence[i] != '.' and sequence[i] != '-':
                matches += 1
        ID = matches/alignment_length


    return

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 HMMER_fasta_summary.py <input_fasta_MSA_file>")
        sys.exit(1)
MSA_file = sys.argv[1]
#output_file = os.path.basename(MSA_file)[0] #work on this
fasta_MSA_summary(MSA_file)
