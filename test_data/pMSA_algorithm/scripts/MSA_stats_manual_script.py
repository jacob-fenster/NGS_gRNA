import os, sys, math, statistics, glob, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import defaultdict
import concurrent, multiprocessing
from concurrent.futures import ProcessPoolExecutor
import subprocess
import time


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a list of descs and a list of seqs.
    """
    descs = []
    seqs = []
    sequence = ""

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()  # Remove leading and trailing whitespaces
            if not line:
                continue  # Skip empty lines
            if line.startswith(">"):  # Description line
                if sequence:
                    seqs.append(sequence)
                    sequence = ""
                descs.append(line[1:])
            else:
                sequence += line

        if sequence:
            seqs.append(sequence)

    return descs, seqs

def query_ID_cov_gap_score(aln, query_aln):
    """ 
    input the alnignment and query alignment in fasta format and output %ID, %Coverage, and %gap score
    %ID = 100 * identical columns/aligned columns, %Coverage = 100 * aligned columns/query length
    %gap = 100 * gaps/query length where query length is number of gap-less columns in query (ignores insertion columns in target)
    """
    aligned, identical, gaps, query_len = 0, 0, 0, 0
    for a, b in zip(aln, query_aln):
        if a != '-' and b != '-':  # Skip gaps 
            aligned += 1
            query_len += 1
            if a == b:  # Count identical positions
                identical += 1
        elif a == '-' and b != '-': #count gaps when present in target column and not in query column
                gaps += 1
                query_len += 1 #query length does not count insertion rows in target
    if aligned == 0 and query_len != 0:
        return 0.0, 100*(aligned/query_len), 100*(gaps/query_len)
    elif aligned == 0 and query_len == 0:
        return 0.0, 0.0, 0.0 
    return 100*(identical/aligned), 100*(aligned/query_len), 100*(gaps/query_len)

def residue_coverage_values(alns, query_aln):
    """
    This one inputs a list of fasta formatted sequence alignments and outputs number of sequences that have residues
    at each position in the query, first aln
    """
    aln_len = len(query_aln)
    residues_cov_x = [i for i in range(aln_len)] 
    residues_cov_y = [0 for i in range(aln_len)]
    for aln in alns:
        for c in range(aln_len):
            if aln[c].isalpha():
                residues_cov_y[c] += 1
    return residues_cov_x, residues_cov_y

def combined_hist_scatter_plt(IDs_perct, covs_perct, residues_cov_x, residues_cov_y, Neff, basename, output_dir):
    tick_fontsize, label_fontsize = 18, 20
    num_values = len(IDs_perct)
    if num_values < 100:
        n_bins = 3 * int(math.sqrt(num_values))
    else: 
        n_bins = int(math.sqrt(num_values))
    # Create a figure and a 1x2 grid of subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    # Add an overall figure title
    title = basename.split('-pMSA')[0]
    fig.suptitle(f'{title} Neff={Neff:.1f} Ntotal={num_values}', fontsize=label_fontsize)
    # Left subplot: Histogram
    axs[0].hist(IDs_perct, bins=n_bins, alpha=0.5, label='%ID', color='b', density=True)
    axs[0].hist(covs_perct, bins=n_bins, alpha=0.6, label='%Coverage', color='#F28C28', density=True)
    axs[0].set_xlabel('%ID or %Coverage', fontsize=label_fontsize)
    axs[0].set_ylabel('Frequency', fontsize=label_fontsize)
    axs[0].set_xlim([0, 100])
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), frameon=False, fontsize=(tick_fontsize-4), framealpha=0)
    # Right subplot: Scatter Plot
    axs[1].plot(residues_cov_x, residues_cov_y, 'k', linewidth=1, label='Residue Coverage')
    axs[1].set_xlabel('Position', fontsize=label_fontsize)
    axs[1].set_ylabel('Sequences', fontsize=label_fontsize)
    axs[1].set_xlim([0, max(residues_cov_x)])
    axs[1].set_ylim([0, max(residues_cov_y)+20])
    #update axis tick fonts
    for tick in axs[0].xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in axs[0].yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in axs[1].xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in axs[1].yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    # Remove top and right spines for left subplot
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    # Increase linewidth and tick size
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)
    axs[0].tick_params(axis='both', width=2, size=6)
    axs[1].spines['left'].set_linewidth(2)
    axs[1].spines['bottom'].set_linewidth(2)
    axs[1].tick_params(axis='both', width=2, size=6)
    # Show the plot
    plt.rcParams.update({'font.size': 24})
    plt.tight_layout() #rect=[0, 0, 1, 1]
    plt.subplots_adjust(top=0.75)
    plt.savefig(f"{output_dir}/{basename}_plt.png")
    plt.close()

def MSA_summary(MSA_list, Nseq_dict, output_dir, fig=False):
    """
    Inputs a list of paths to MSAs or pMSAs in fasta format, calculates query specific stats, Neff sequences, and adds to a row in a growing dataframe
    First line of MSA/pMSA must be query alignment.
    with the MSA speicifc stats it generates a histogram of the %ID and %Coverage distribution with summary stats.  
    """
    summary = pd.DataFrame(columns=['#_alns', 'Nseq_80', 'Neff_80', 'mean_ID', 'stdev_ID', 'mean_cov', 'stdev_cov'])
    for MSA in MSA_list:
        index = os.path.basename(MSA).split('.')[0] #name the summary row by the MSA/pMSA basename
        descs, alns = read_fasta(MSA)
        IDs, covs = [], []
        query_aln = alns[0] #first entry must be query sequence
        for i in range(len(alns)): #find vs query summary stats
            ID_perct, cov_perct, gap_perct = query_ID_cov_gap_score(alns[i], query_aln)
            IDs.append(ID_perct)
            covs.append(cov_perct)
        residues_cov_x, residues_cov_y = residue_coverage_values(alns, query_aln)
        Nseq, Neff = Nseq_dict[index], Nseq_dict[index]/math.sqrt(len(query_aln))
        summary.loc[index, '#_alns'], summary.loc[index, 'Nseq_80'], summary.loc[index, 'Neff_80'] = len(alns), Nseq, Neff
        summary.loc[index, 'mean_ID'], summary.loc[index, 'mean_cov'] = statistics.mean(IDs), statistics.mean(covs)
        if len(IDs) < 2:
            summary.loc[index, 'stdev_ID'], summary.loc[index, 'stdev_cov'] = 0, 0
        else:
            summary.loc[index, 'stdev_ID'], summary.loc[index, 'stdev_cov'] = statistics.stdev(IDs), statistics.stdev(covs)
        if fig:
            combined_hist_scatter_plt(IDs, covs, residues_cov_x, residues_cov_y, Neff, index, output_dir)
    return summary

def MSA_summary_AFNeff(MSA_list, MSA_clust_list, output_dir, fig=False):
    """
    Inputs a list of paths to MSAs or pMSAs in fasta format, with corresponding hhfilter clustered MSAs list in fasta format,
    calculates query specific stats, Neff sequences, and adds to a row in a growing dataframe
    First line of MSA/pMSA must be query alignment.
    with the MSA speicifc stats it generates a histogram of the %ID and %Coverage distribution with summary stats.  
    """
    summary = pd.DataFrame(columns=['#_alns', 'Nseq_80', 'Neff_80', 'AF_Neff_80', 'mean_ID', 'stdev_ID', 'mean_cov', 'stdev_cov'])
    for MSA, MSA_clust in zip(MSA_list, MSA_clust_list):
        index = os.path.basename(MSA).split('.')[0] #name the summary row by the MSA/pMSA basename
        descs, alns = read_fasta(MSA)
        descs_clust, alns_clust = read_fasta(MSA_clust)
        IDs, covs = [], []
        query_aln = alns[0] #first entry must be query sequence
        for i in range(len(alns)): #find vs query summary stats
            ID_perct, cov_perct, gap_perct = query_ID_cov_gap_score(alns[i], query_aln)
            IDs.append(ID_perct)
            covs.append(cov_perct)
        residues_cov_x, residues_cov_y = residue_coverage_values(alns, query_aln)
        residues_cov_x_clust, residues_cov_y_clust = residue_coverage_values(alns_clust, query_aln)
        Nseq, Neff = len(descs_clust), len(descs_clust)/math.sqrt(len(query_aln))
        AF_Neff = statistics.median(residues_cov_y_clust)
        summary.loc[index, '#_alns'], summary.loc[index, 'Nseq_80'], summary.loc[index, 'Neff_80'] = len(alns), Nseq, Neff
        summary.loc[index, "AF_Neff_80"] = AF_Neff
        summary.loc[index, 'mean_ID'], summary.loc[index, 'mean_cov'] = statistics.mean(IDs), statistics.mean(covs)
        if len(IDs) < 2:
            summary.loc[index, 'stdev_ID'], summary.loc[index, 'stdev_cov'] = 0, 0
        else:
            summary.loc[index, 'stdev_ID'], summary.loc[index, 'stdev_cov'] = statistics.stdev(IDs), statistics.stdev(covs)
        if fig:
            combined_hist_scatter_plt(IDs, covs, residues_cov_x, residues_cov_y, Neff, index, output_dir)
    return summary

def four_plot_figure(num_alns, avg_IDs, avg_covs, N_80, Neff, exp_tag, output_dir):
    tick_fontsize, label_fontsize = 22, 24
    num_values = len(avg_IDs)
    if num_values < 100:
        n_bins = 3 * int(math.sqrt(num_values))
    else: 
        n_bins = int(math.sqrt(num_values))
        
    # Create a figure and a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))
    
    # Add an overall figure title
    title = exp_tag
    fig.suptitle(f'{title}', fontsize=label_fontsize+2)
    
    # Top-left subplot: Histogram of #_alns
    axs[0, 0].hist(num_alns, bins=n_bins, color='g', density=True)
    axs[0, 0].set_xlim([0, None])
    axs[0, 0].set_xlabel('#_alns', fontsize=label_fontsize)
    axs[0, 0].set_ylabel('Frequency', fontsize=label_fontsize)
    current_ticks = axs[0, 0].get_xticks()
    new_ticks = current_ticks[::2]  # take every second tick
    axs[0, 0].set_xticks(new_ticks)

    # Top-right subplot: Histogram of %ID and %Coverage
    axs[0, 1].hist(avg_IDs, bins=n_bins, alpha=0.5, label='average %ID', color='b', density=True)
    axs[0, 1].hist(avg_covs, bins=n_bins, alpha=0.6, label='average %Coverage', color='#F28C28', density=True)
    axs[0, 1].set_xlabel('avg %ID or avg %cov', fontsize=label_fontsize)
    axs[0, 1].set_xlim([0, 100])
    axs[0, 1].legend(loc='upper center', fontsize=(tick_fontsize-4), frameon=False)
    
    # Bottom-left subplot: Histogram of N_80
    axs[1, 0].hist(N_80, bins=n_bins, color='r', density=True)
    axs[1, 0].set_xlim([0, None])
    axs[1, 0].set_xlabel('N_80', fontsize=label_fontsize)
    axs[1, 0].set_ylabel('Frequency', fontsize=label_fontsize)
    
    # Bottom-right subplot: Histogram of Neff
    axs[1, 1].hist(Neff, bins=n_bins, color='m', density=True)
    axs[1, 1].set_xlim([0, None])
    axs[1, 1].set_xlabel('Neff', fontsize=label_fontsize)
    axs[1, 1].set_ylabel('Frequency', fontsize=label_fontsize)
    
    # Update axis tick fonts and style for all subplots
    for i in range(2):
        for j in range(2):
            for tick in axs[i, j].xaxis.get_major_ticks():
                tick.label1.set_fontsize(tick_fontsize)
            for tick in axs[i, j].yaxis.get_major_ticks():
                tick.label1.set_fontsize(tick_fontsize)
            axs[i, j].spines['top'].set_visible(False)
            axs[i, j].spines['right'].set_visible(False)
            axs[i, j].spines['left'].set_linewidth(2)
            axs[i, j].spines['bottom'].set_linewidth(2)
            axs[i, j].tick_params(axis='both', width=2, size=6)

    #plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(f"{output_dir}/{exp_tag}_runsummaryplt.png")
    plt.close()

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

# input

File = False
a3m_file = "/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/FastaBank_clustalo/pMSAs/F778R-AA778__F334L-AA334-pMSA.a3m"

directory = True
input_a3m_dir = '/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/clustered_at_99/pMSAs_clust99'

temp_dir = '/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/clustered_at_99/temp'
output_dir = '/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/clustered_at_99/stats'
exp_tag = 'ASFV-ASFV_Clust99_AF_Neff80'

a3m_files = glob.glob(f"{input_a3m_dir}/*.a3m")
hhfilter = '/Users/jacobfenster/Documents/ASFV_Postdoc/hhsuite/bin/hhfilter'
reformat = '/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/scripts/reformat.pl'
Nseq_dict = {}
MSA_list, clust_MSA_list = [], []
no_file_error, no_alns = 0, 0
test_parallel = True
if __name__ == '__main__':
    if directory == True:
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        clust_commands = []
        reformat_commands = []
        reformat_clust_commands = []
        for file in a3m_files:
            basename = os.path.basename(file).split('.')[0]
            clust_commands.append(f"{hhfilter} -v 0 -id 80 -i {file} -o {temp_dir}/{basename}_clust80.a3m")
            reformat_commands.append(f"{reformat} -v 0 a3m fas {file} {temp_dir}/{basename}.fasta")
            reformat_clust_commands.append(f"{reformat} -v 0 a3m fas {temp_dir}/{basename}_clust80.a3m {temp_dir}/{basename}_clust80.fasta")
        timei = time.time()
        results_clust = run_pooled_commands(clust_commands)
        results_reformat = run_pooled_commands(reformat_commands)
        clust_reformat = run_pooled_commands(reformat_clust_commands)
        timef = time.time()
        print(f'there has elapsed {timef-timei} seconds')
        for file in a3m_files:
            basename = os.path.basename(file).split('.')[0] 
            MSA_list.append(f"{temp_dir}/{basename}.fasta")
            clust_MSA_list.append(f"{temp_dir}/{basename}_clust80.fasta")
        summary = MSA_summary_AFNeff(MSA_list, clust_MSA_list, output_dir, fig=False)
        summary.to_csv(f"{output_dir}/{exp_tag}_summary.csv")
        os.system(f"rm -r {temp_dir}")
    elif False:
        os.mkdir(temp_dir)
        os.mkdir(output_dir, exist_ok=True)
        for file in a3m_files:
            basename = os.path.basename(file).split('.')[0]
            os.system(f"{hhfilter} -v 0 -id 90 -i {file} -o {temp_dir}/{basename}_clust80.a3m")
            os.system(f"{reformat} -v 0 a3m fas {file} {temp_dir}/{basename}.fasta")
            try:
                descs, seqs = read_fasta(f"{temp_dir}/{basename}_clust80.a3m")
                if len(descs) == 0:
                    no_alns += 1
                    continue
                Nseq_dict[basename] = len(descs)
                MSA_list.append(f"{temp_dir}/{basename}.fasta")
            except FileNotFoundError:
                no_file_error += 1 
                continue
        print(f"there are {no_file_error} no file errors for this run. There are {no_alns} pMSAs with zero alns.")
        summary = MSA_summary(MSA_list, Nseq_dict, output_dir)
        summary.to_csv(f"{output_dir}/{exp_tag}_summary.csv")
        num_alns = summary['#_alns'].tolist()
        avg_IDs = summary['mean_ID'].tolist()
        avg_covs = summary['mean_cov'].tolist()
        N_80 = summary['Nseq_80'].tolist()
        Neff = summary['Neff_80'].tolist()
        #four_plot_figure(num_alns, avg_IDs, avg_covs, N_80, Neff, exp_tag, output_dir)  
        os.system(f"rm -r {temp_dir}")
        pdb.set_trace()


    elif File == True:
        os.mkdir(temp_dir)
        basename = os.path.basename(a3m_file).split('.')[0]
        os.system(f"{hhfilter} -v 0 -id 80 -i {a3m_file} -o {temp_dir}/{basename}_clust80.a3m")
        os.system(f"{reformat} a3m fas {a3m_file} {temp_dir}/{basename}.fasta")
        descs, seqs = read_fasta(f"{temp_dir}/{basename}_clust80.a3m")
        Nseq_dict[basename] = len(descs)
        MSA_list.append(f"{temp_dir}/{basename}.fasta")
        summary = MSA_summary(MSA_list, Nseq_dict, output_dir, fig=True)
        summary.to_csv(f"{output_dir}/{exp_tag}_summary.csv")
        #os.system(f"rm -r {temp_dir}")
#there has elapsed 54.47777199745178 seconds