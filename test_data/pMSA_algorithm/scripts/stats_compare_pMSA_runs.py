import os
import sys 
import pandas as pd
import numpy as np
import math
import statistics
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import defaultdict
import pdb

def four_plot_figure(num_alns, avg_IDs, avg_covs, N_80, Neff, basename):
    tick_fontsize, label_fontsize = 22, 24
    num_values = len(avg_IDs)
    if num_values < 100:
        n_bins = 3 * int(math.sqrt(num_values))
    else: 
        n_bins = int(math.sqrt(num_values))
        
    # Create a figure and a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))
    
    # Add an overall figure title
    title_basename = os.path.basename(basename).split('-pMSA')[0]
    fig.suptitle(f'{title_basename}', fontsize=label_fontsize+2)
    
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
    plt.savefig(f"{basename}_runsummaryplt.png")
    plt.close()



run, output = sys.argv[1], sys.argv[2]

#output summary histogram plots of #_aln, N_80, Neff, mean_ID, mean_cov
# #_aln single, N_80 & Neff double, mean_ID & mean_cov double single plot 
summary = pd.read_csv(run)
# Dummy data for testing
num_alns = summary['#_alns'].tolist()
avg_IDs = summary['mean_ID'].tolist()
avg_covs = summary['mean_cov'].tolist()
N_80 = summary['Nseq_80'].tolist()
Neff = summary['Neff_80'].tolist()


# Test the function
four_plot_figure(num_alns, avg_IDs, avg_covs, N_80, Neff, output)






