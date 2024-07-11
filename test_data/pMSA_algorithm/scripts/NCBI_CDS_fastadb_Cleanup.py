import sys
import os
from Bio import SeqIO
import re

def filter_fasta_summary(input_file_path):
    output_file_path = f"{input_file_path.split('.')[0]}_cleanup.fasta"
    removed_file_path = f"{input_file_path.split('.')[0]}_cleanup_dropped_CDS.fasta"
    output_summary_path = f"{input_file_path.split('.')[0]}_cleanup_summary.txt"
    removed_summary_path = f"{input_file_path.split('.')[0]}_cleanup_dropped_CDS_summary.txt"
    with open(input_file_path, "r") as input_file, open(output_file_path, "w") as output_file, open(removed_file_path, "w") as dropped_file:
        unique_genomes_out, unique_proteins_out, unique_genomes_drop, unique_proteins_drop = set(), set(), set(), set()
        for record in SeqIO.parse(input_file, "fasta"):
            header = record.description
            genome_id = header.split('|')[1].split('_prot_')[0].split('.')[0]
            protein_id = header.split('|')[1].split('_prot_')[1].split('.')[0]
            if all(char not in record.seq for char in ['X', '-', 'J', 'B', 'Z']):
                SeqIO.write(record, output_file, "fasta")
                unique_genomes_out.add(genome_id)
                unique_proteins_out.add(protein_id)
            else: 
                SeqIO.write(record, dropped_file, "fasta")
                unique_genomes_drop.add(genome_id)
                unique_proteins_drop.add(protein_id)
    #output summary files
    with open(output_summary_path, "w") as output_summary, open(removed_summary_path, "w") as removed_summary:
        output_summary.write(f"# Original_genomes={len(unique_genomes_out)+len(unique_genomes_drop)}\tOriginal_proteins={len(unique_proteins_out)+len(unique_proteins_drop)}\n")
        output_summary.write(f"# Genomes_after_filtering={len(unique_genomes_out)}\tProteins_after_filtering={len(unique_proteins_out)}\n# Genome_IDs_included\n")
        removed_summary.write(f"# Removed_genomes={len(unique_genomes_drop)}\tRemoved_proteins={len(unique_proteins_drop)}\n# Genome_IDs_removed\n")
        for genome in unique_genomes_out:
            output_summary.write(f"{genome}\n")
        for genome in unique_genomes_drop:
            removed_summary.write(f"{genome}\n")

input_file_path = sys.argv[1]
filter_fasta_summary(input_file_path)