import ipdb
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def plot_replicates(df, condition, axis_limits=None, ticks=None):
    x = df[f"{condition}_R1"]
    y = df[f"{condition}_R2"]

    # Scatter plot
    plt.figure(figsize=(3, 3))
    sns.scatterplot(x=x, y=y, s=1, color='black', edgecolor=None)
    plt.xlabel("Rep 1")
    plt.ylabel("Rep 2")
    plt.title(f"{condition}")
    # Set axis limits if provided
    if axis_limits is not None:
        plt.xlim(axis_limits)
        plt.ylim(axis_limits)
    if ticks is not None:
        plt.xticks(ticks)
        plt.yticks(ticks)

    # Linear regression
    reg = LinearRegression().fit(x.values.reshape(-1, 1), y.values.reshape(-1, 1))
    y_pred = reg.predict(x.values.reshape(-1, 1))
    r2 = r2_score(y, y_pred)
    
    # Plot the regression line
    plt.plot(x, y_pred, color='red', label=f"Linear fit (R² = {r2:.2f})")
    plt.legend()
    plt.show()

def plot_histogram(df, column_name, xlim=None):
    """
    Plots a histogram of a specified column from a given DataFrame.
    
    Parameters:
    df (pandas.DataFrame): The DataFrame containing the data.
    column_name (str): The name of the column to plot.
    """
    # Create the histogram
    fig, ax = plt.subplots(figsize=(3, 3))

    ax.hist(df[column_name], bins=100, color='blue', edgecolor='black')

    # Set the x and y labels
    ax.set_xlabel('gRNA reads')
    ax.set_ylabel('#gRNAs')

    # Reduce the number of ticks on x and y axis
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    if xlim is not None:
        plt.xlim(xlim)
    # Show the plot
    plt.tight_layout()
    plt.show()

df = pd.read_csv('/Users/jfenster/Documents/NGS_gRNA/output/20231126_exact_counts/20231126_Carbon_exact_MAGeCK.txt', 
                 sep='\t')
# Normalize the counts
plot_histogram(df, '50pCA_t1_R1', xlim=(0,2000))
count_columns = df.columns[2:]
for column in count_columns:
    df[column] = df[column] / df[column].sum()
df.fillna(0, inplace=True)
# Generate scatter plots for each condition
conditions = ['LB_t0', 'LB_t1', 'LB_t2', 'gluc_t1', 'gluc_t2', 'ace_t1', 'ace_t2', '10pCA_t1', '10pCA_t2', '50pCA_t1', '50pCA_t2']
breakpoint()
plot_replicates(df, '50pCA_t1', axis_limits=(0, .0001), ticks=[0,.0001])
breakpoint()
for condition in conditions:
    plot_replicates(df, condition, axis_limits=.0001)
    breakpoint()