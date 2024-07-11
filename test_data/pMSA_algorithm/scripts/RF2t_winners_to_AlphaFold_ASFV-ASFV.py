import os
import pandas as pd
import sys

def read_fasta_dict(fasta_file):
    entries = {}
    sequence = ""
    header = None
    with open(fasta_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries[header] = sequence
                header = line[1:]
                sequence = ""
            else:
                sequence += line
        if header:
            entries[header] = sequence
    return entries


#Usage: $python3 RF2t_winners_to_AlphaFold_ASFV-ASFV.py <path/to/ASFV genome.fasta> <path/to/output/dir> </path/to/RF2t++.csv data file> <RF2t++ Threshold cuttoff>

ASFV_genome_file, output_dir, RF2t_csv_path, RF2t_plusplus_tr = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
genome = os.path.basename(ASFV_genome_file).split('.')[0]
entries = read_fasta_dict(ASFV_genome_file)
RF2t_df = pd.read_csv(RF2t_csv_path)
winners_df = RF2t_df[RF2t_df['RF2t++ scores'] > RF2t_plusplus_tr]
winners = winners_df.sort_values(by='RF2t++ scores', ascending=False)['PPI name'].tolist()

if '-AA' in winners[0]: #switch to catch the formatting Protein1-AA###__Protein2-AA###
    index = 0
    for pair in winners:
        protein1, protein2 = pair.split('-AA')[0], pair.split('__')[1].split('-AA')[0]
        if "NoQuery" in protein1 or "NoQuery" in protein2: #skip proteins that are not present in ASFV-G 2007
            index +=1
            continue
        if "ACD" in protein1:
            protein1 = "ASFV_G_"+protein1
        if "ACD" in protein2:
            protein2 = "ASFV_G_"+protein2
        if protein2 == 'MGF_110-10-L-MGF110-14L_fusion':
            protein2 = "MGF_110-10-L_-_MGF110-14L_fusion"
        with open(f"{output_dir}/{str(index).zfill(3)}_ASFVG_{protein1}__{protein2}.fa", "w") as outfile:
            outfile.write(f">{protein1}|{genome}\n")
            for i in range(0, len(entries[protein1]), 60):
                outfile.write(f"{entries[protein1][i:i+60]}\n")
            outfile.write(f"\n>{protein2}|{genome}\n")
            for i in range(0, len(entries[protein2]), 60):
                outfile.write(f"{entries[protein2][i:i+60]}\n")
        index += 1
        
else: 
    print("This script does not support that file format")
             

