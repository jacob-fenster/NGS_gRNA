import os
import sys
from Bio import SeqIO

def split_proteome(input_file_path):
    for record in SeqIO.parse(input_file_path, 'fasta'):
        #this just puts the fasta entry description as the file name for output. Changes the protein name to query for pMSA
        #outputs to working directory as well for nextflow usage
        output_file = f"{record.description}.fasta"
        with open(output_file, 'w') as f:
            f.write(f">query\n{str(record.seq)}\n")

input_file_path = sys.argv[1]
split_proteome(input_file_path)