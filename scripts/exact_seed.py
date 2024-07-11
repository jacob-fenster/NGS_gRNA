import sys, os, pdb
from Bio.Seq import Seq
import pandas as pd
from collections import Counter


def extract_exact_seed(merged_fastq_filepath, up_motif, down_motif):
    seeds, not_exact, not_exact_len = [], [], []
    n, reads = 0, 0
    with open(merged_fastq_filepath, 'r') as f:
        for line in f:
            if (n - 1) % 4 == 0: #read each second line of each input
                reads += 1
                seq = Seq(line.strip())
                seq_rc = seq.reverse_complement()
                up_index, down_index = seq.find(up_motif), seq.find(down_motif)
                up_index_rc, down_index_rc = seq_rc.find(up_motif), seq_rc.find(down_motif)
                if up_index != -1 and down_index != -1: # might need to put case where there is more than one index
                    seed = str(seq[up_index+35:down_index])
                    if len(seed) == 20:
                        seeds.append(seed)
                elif up_index_rc != -1 and down_index_rc != -1:
                    seed = str(seq_rc[up_index_rc+35:down_index_rc])
                    if len(seed) == 20:
                        seeds.append(seed)
                else:
                    not_exact.append(str(seq))
                    not_exact_len.append(len(str(seq)))
            n += 1
    return reads, seeds, not_exact, not_exact_len

def exact_read_counts(seed_reads, seed_db_path):
    pdb.set_trace()
    seed_df = pd.read_csv(seed_db_path, index_col=0)
    seed_db, gene_db, name_db = seed_df['seed'].tolist(),seed_df['gene'].tolist(), seed_df.index.tolist()
    seed_dict = dict(zip(seed_db, name_db))
    seed_read_counts = Counter(seed_reads)
    reads = seed_df.copy()
    incorrect = 0
    for seed, count in seed_read_counts.items():
        try:
            reads.loc[seed_dict[seed], 'reads'] = count
        except KeyError:
            incorrect += 1
    return reads, incorrect


merged_fastq_filepath, seed_db_path, output_filename = sys.argv[1], sys.argv[2], sys.argv[3]
#these are specific to the 78K P. putida KT2440 gNRA libraries
#after trimming the up and down ovlp sequences below are the longest segment that contains the overlapping region 5' and 3'
#the primers were phased, so the start varies by 4 bp for correct sequences. N20 spacer/barcode is in between these
up_ovlp_motif = Seq('TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC')
dwn_ovlp_motif = Seq('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAA')
reads, seeds, not_exact, not_exact_len = extract_exact_seed(merged_fastq_filepath, up_ovlp_motif, dwn_ovlp_motif)
hits, fails = len(seeds), len(not_exact)
print(f"Reads: {reads} Exact hit rate: {int(100*hits/reads)}% for file {os.path.basename(merged_fastq_filepath)}")
pdb.set_trace()
reads, incorrect = exact_read_counts(seeds, seed_db_path)
reads.to_csv(output_filename)
pdb.set_trace()