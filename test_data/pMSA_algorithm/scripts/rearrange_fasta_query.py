import sys
import os
import textwrap  # Importing the textwrap library for easy string formatting

# Function to format sequence with 60 characters per line
def format_sequence(sequence):
    return "\n".join(textwrap.wrap(sequence, 60))

# This rearranges the input fasta placing the entry that matches the query on top
# It also outputs a single entry query fasta file for downstream search/MSA
def rearrange_fasta(input_file, query_string, output_file, query_file):
    entries = {}
    sequence = ""
    header = None

    # Read the FASTA file and populate the 'entries' dictionary
    with open(input_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries[header] = sequence
                header = line
                sequence = ""
            else:
                sequence += line
        if header:
            entries[header] = sequence

    # Search for the query and place it at the top
    top_entry = None
    for header, sequence in entries.items():
        if query_string in header:
            top_entry = (header, sequence)
            del entries[header]
            break

    if top_entry is None: #write to ouput anyway with augmented name if query not found
        # Write to the output file
        output_file = 'NoQuery'+output_file
        with open(output_file, "w") as outfile:
            for header, sequence in entries.items():
                outfile.write(f"{header}\n{format_sequence(sequence)}\n")
        
    else:
        # Write to the output file
        with open(output_file, "w") as outfile:
            outfile.write(f"{top_entry[0]}\n{format_sequence(top_entry[1])}\n")
            for header, sequence in entries.items():
                outfile.write(f"{header}\n{format_sequence(sequence)}\n")

        # Write the query fasta output file
        with open(query_file, "w") as outfile:
            outfile.write(f"{top_entry[0]}\n{format_sequence(top_entry[1])}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python rearrange_fasta.py <input_fasta_file> <query_string> <output_fasta_filename> <output_fasta_query_filename>")
        sys.exit(1)


    input_file = sys.argv[1]
    query_string = sys.argv[2]
    output_file = sys.argv[3]
    query_file = sys.argv[4]

    rearrange_fasta(input_file, query_string, output_file, query_file)
