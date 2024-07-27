import os, ipdb
import pandas as pd
import math
import statistics 
import seaborn as sns
import matplotlib.pyplot as plt

def plot_with_error_bars(x, y, error):
    """
    Plots a graph with error bars using seaborn and matplotlib.

    Parameters:
    x (list or array): x values
    y (list or array): y values
    error (list or array): error values for y
    """
    # Create a DataFrame from your data
    data = pd.DataFrame({'x': x, 'y': y, 'error': error})

    # Set font size
    plt.rcParams.update({'font.size': 12})

    # Create the plot with specified figure size
    fig, ax = plt.subplots(figsize=(3.5, 2.5))

    # Plot the points with error bars
    ax.errorbar(data['x'], data['y'], yerr=data['error'], fmt='o', color='black', ecolor='black', capsize=5, markersize=5)

    # Customize the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add x and y axis lines
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')

    # Set x-axis limit
    ax.set_xlim(left=-25)

    # Display the plot
    plt.show()
data_dir = "/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/"
pca_data = pd.read_csv(f"{data_dir}50pCA_t2-vs-50pCA_t1.sgrna_summary.txt", 
                       sep='\t', index_col=0)
mreb_data = pca_data[pca_data['Gene'] == 'mreB']
# calculate relative cds location
locations = []
for index in mreb_data.index:
    start = int(index.split(':')[1].split('-')[0])
    locations.append(start)
zero_index = min(locations)
variance, distance, pct_mu = [], [], []

for index in mreb_data.index:
    start = int(index.split(':')[1].split('-')[0])
    location =  start - zero_index
    distance.append(location)
    lfc = mreb_data.loc[index, 'LFC']
    pct_mu.append((lfc/10 ** .5))
    # pull out counts 
    c1, c2 = float(mreb_data.loc[index, 'control_count'].split('/')[0]), float(mreb_data.loc[index, 'control_count'].split('/')[1])
    t1, t2 = float(mreb_data.loc[index, 'treatment_count'].split('/')[0]), float(mreb_data.loc[index, 'treatment_count'].split('/')[1])
    mean_c, mean_t = statistics.mean([c1, c2]), statistics.mean([t1, t2])
    std_c, std_t = statistics.stdev([c1, c2]), statistics.stdev([t1, t2])
    error_div = lfc * math.sqrt((std_c/mean_c) ** 2 + (std_t/mean_t) ** 2)
    error = (error_div / (lfc * math.log(2))) / 10 ** .5
    variance.append(error)
distance = [i + 25 for i in distance]
plot_with_error_bars(distance, pct_mu, variance)
breakpoint()



