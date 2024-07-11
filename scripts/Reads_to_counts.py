import os, sys, glob, re
import pandas as pd
from Bio.Seq import Seq
from collections import Counter
import concurrent, multiprocessing
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import subprocess

def extract_exact_seed(merged_fastq_filepath, up_motif, down_motif):
    """
    extracts only 20 bp seed regions that have EXACTLY the up and down motif
    returns the number of reads, reads. a list of seeds for each read, seed. 
    a list of reads that were not exact for troubleshooting, not_exact, and the length of each not exact read, not_exact_len
    Also returns the basename to pair file-results for parallel processing. Must have "_merged" tag
    """
    sample_name = os.path.basename(merged_fastq_filepath).split('_merged')[0] #specific to _merged tag on file
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
                    else: 
                        not_exact.append(str(seq))
                elif up_index_rc != -1 and down_index_rc != -1:
                    seed = str(seq_rc[up_index_rc+35:down_index_rc])
                    if len(seed) == 20:
                        seeds.append(seed)
                    else: 
                        not_exact.append(str(seq))
                else:
                    not_exact.append(str(seq))
            n += 1
    return sample_name, seeds, reads, len(not_exact)

def exact_read_counts(seeds, seed_df, sample_name):
    """
    This inputs seeds from extract_exact_seed and the path to the sgRNA databse 
    outputs a read dataframe and a count of incorrect seeds that are not in the dict
    """
    seed_db, gene_db, name_db = seed_df['seed'].tolist(), seed_df['gene'].tolist(), seed_df.index.tolist()
    seed_dict = dict(zip(seed_db, name_db))
    seed_read_counts = Counter(seeds)
    counts_df = seed_df.copy()
    counts_df[sample_name] = 0 #initialize reads to zero
    incorrect = 0
    for seed, count in seed_read_counts.items():
        try:
            counts_df.loc[seed_dict[seed], sample_name] = count
        except KeyError:
            incorrect += count
    return counts_df, incorrect

def extract_counts_parallel(merged_fastq_filepaths, up_motif, down_motif, seed_df):
    """
    This inputs the list of merged fastq file paths, the up and down motif, and the starting seed_df (database) and calls the extract exact read in parallel
    To simplify summary stats, the exact_read_counts is called in series and these are added ot the counts_summary and read_stats df
    """
    counts_summary = seed_df.copy()
    read_stats = pd.DataFrame(columns=['total_merged_reads', 'exact_reads', '%_exact_reads', 'not_exact_reads', '%_not_exact_reads'])
    num_cores = multiprocessing.cpu_count() - 1
    samples = len(merged_fastq_filepaths)
    if num_cores < 1:
        num_cores = 1
    num_workers = min(len(merged_fastq_filepaths), num_cores)
    up_motifs = repeat(up_motif, samples) #convert constants to iterables
    down_motifs = repeat(down_motif, samples)   
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        seed_results = executor.map(extract_exact_seed, merged_fastq_filepaths, up_motifs, down_motifs)
    # Pull out results for parallel counts
    seed_sample_list, seeds_lists, total_reads, not_exact_reads = [], [], [], []
    for result in seed_results:
        sample_name, seeds, reads, not_exact = result
        seed_sample_list.append(sample_name), seeds_lists.append(seeds), total_reads.append(reads), not_exact_reads.append(not_exact)
    seed_df_repeat = repeat(seed_df, samples)
    for seeds, seed_df, sample, total_read, not_exact_read in zip(seeds_lists, seed_df_repeat, seed_sample_list, total_reads, not_exact_reads):
        counts_df, incorrect = exact_read_counts(seeds, seed_df, sample)
        counts_summary[sample] = counts_df[sample]
        total_incorrect, exact_reads = incorrect + not_exact_read, counts_df[sample].sum()
        read_stats.loc[sample] = {'total_merged_reads':total_read, 'exact_reads':exact_reads, '%_exact_reads':100*exact_reads/total_read, 'not_exact_reads':total_incorrect, '%_not_exact_reads':100*total_incorrect/total_read}
    return counts_summary, read_stats

def find_file_pairs(directory):
    """
    This is specific to finding pairs of files named *_R[1,2]_*.fastq
    """
    # Create a dictionary to store file pairs
    file_pairs = {}

    # Pattern to match the files (assumes the format is like R1-1_R1_001.fastq)
    pattern = re.compile(r'(.+)_R[12]_.*\.fastq')

    # List all .fastq files in the directory
    for file in glob.glob(os.path.join(directory, '*_R[12]_001.fastq')):
        # Extract the base name without R1 or R2
        match = pattern.match(os.path.basename(file))
        if match:
            base_name = match.group(1)

            # Add the file to the corresponding pair in the dictionary
            if base_name in file_pairs:
                file_pairs[base_name].append(file)
            else:
                file_pairs[base_name] = [file]

    # Filter out any unmatched pairs and return
    return {k: v for k, v in file_pairs.items() if len(v) == 2}

def run_command(command):
    result = subprocess.run(command, shell=True)
    return result.returncode

def run_pooled_commands(commands):
    """
    Inputs a list of system commands and runs them on parallel CPU cores. 
    Determines the number of CPU cores to use by 1 minus the available cores on the system. 
    Defaults to 1 CPU core in case resources are limited upon execution.
    """
    num_cores = multiprocessing.cpu_count() - 1
    if num_cores < 1:
        num_cores = 1
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        results = executor.map(run_command, commands)
    return list(results)

if __name__ == '__main__':
    reads_folder = '/Volumes/Extreme_SSD/78K_CRISPRi/Expanded'
    temp_dir = '/Volumes/Extreme_SSD/78K_CRISPRi/Expanded/work'
    output_dir = '/Users/jacobfenster/Documents/NGS_gRNA/output/20231126'
    seed_db_path = "/Users/jacobfenster/Documents/NGS_gRNA/databases/Lib3_dCas9_78K_seeds.csv"
    exp_tag = '20231126_carbon_sources'

    up_motif = Seq('TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC')
    down_motif = Seq('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAA')
    bbmerge = '/Users/jacobfenster/Documents/NGS_gRNA/docker_containers/bbmerge/software/bbmap/bbmerge.sh'

    seed_df = pd.read_csv(seed_db_path, index_col=0)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    pairs = find_file_pairs(reads_folder)
    merge_commands = []
    for pair in pairs.values():
        sampleID = os.path.basename(pair[0]).split('_R')[0]
        merge_commands.append(f"{bbmerge} in={pair[0]} in2={pair[1]} out={temp_dir}/{sampleID}_merged.fastq ihist=${output_dir}/{sampleID}_merge_hist.txt trimnonoverlapping=t")
    #run_pooled_commands(merge_commands)
    breakpoint()
    merged_fastq = glob.glob(f"{temp_dir}/*_merged.fastq")
    counts_summary, read_stats = extract_counts_parallel(merged_fastq, up_motif, down_motif, seed_df)
    counts_summary.to_csv(f"{output_dir}/{exp_tag}-counts.tsv", sep='\t')
    read_stats.to_csv(f"{output_dir}/{exp_tag}-readstats.tsv", sep='\t')
    # sampleID = 'R1-6'
    # f"{bbmerge} in={pairs['R1-6'][0]} in2={pairs['R1-6'][1]} out={temp_dir}/{sampleID}_merged.fastq ihist=${output_dir}/{sampleID}_merge_hist.txt trimnonoverlapping=t"

