import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pdb

def plot_histograms_NormPercent(df, col, file_basename):
    """
    Plots histograms for specififed column
    - df (DataFrame): The input DataFrame containing the data
    - column_names (list): List of column names for which to plot histograms
    """
    #convert columns of interest to floats
    df[col] = df[col].astype(float)
    # Calculate statistics
    mean_value, median_value, stdev, num_values = df[col].mean(), df[col].median(), df[col].std(), df[col].count()
    if num_values < 100:
        n_bins = 3 * int(np.sqrt(len(df[col])))
    else: 
        n_bins = int(np.sqrt(len(df[col])))
    plt.figure(figsize=(6, 6))
    plt.hist(df[col], bins=n_bins, edgecolor='black', alpha=0.7, density=True)
    plt.xlim(0,100)
    plt.rcParams.update({'font.size': 20})
    plt.annotate(f"Mean: {mean_value:.1f}\nMedian: {median_value:.1f}\nStd Dev: {stdev:.1f}\nn: {num_values}\nn_bins: {n_bins}", xy=(0.6, 0.8), xycoords='axes fraction', fontsize=12)
    plt.title(file_basename+"_"+col)
    plt.xlabel(col)
    plt.ylabel('Frequency')
    # Save the plot
    plt.tight_layout()
    plt.savefig(file_basename+'_'+col+'_hist.png')
    plt.close()
    return


#Usage python3 $projectDir/scripts/hist_normpercent.py "${MSA_csv}" "${columns}"
file_path = sys.argv[1]
column = sys.argv[2] #names of columns to generate histograms of
df = pd.read_csv(file_path)
if '10239' in file_path: #check to see if orthoDB ASFV level. specific to my nonsense naming scheme
    file_basename = os.path.basename(file_path).split('__')[0].split('_')[1] + os.path.basename(file_path).split('.')[0].split('10239')[1].split('_')[0]
else: #catch all for other file types
    file_basename = os.path.basename(file_path).split('.')[0].split('_rearranged_clustalo_filtered')[0]

if column not in df.columns:
    print(f"Column '{col}' not found in DataFrame. Skipping...")
else: 
    plot_histograms_NormPercent(df, column, file_basename)