import csv
import sys
import os

# Function to convert CSV to FASTA
def csv_to_fasta(csv_file, fasta_file):
    with open(csv_file, 'r') as csv_f:
        reader = csv.reader(csv_f)
        with open(fasta_file, 'w') as fasta_f:
            for row in reader:
                description = row[0]
                sequence = row[1]
                fasta_f.write(f">{description}\n{sequence}\n")

# Input CSV file
csv_file = sys.argv[1]

# Output FASTA file
file_name, _ = os.path.splitext(csv_file)
fasta_file = file_name + '.fasta'

# Convert CSV to FASTA
csv_to_fasta(csv_file, fasta_file)
print(f"Converted {csv_file} to {fasta_file}")

#running $ python3 scripts/csv_to_fasta.py databases/PPutida_gRNA_Libraries/Ordered_filtered_databases/Lib3_dCas9_filtered_ordered.csv