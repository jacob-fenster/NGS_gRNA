import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import ipdb

# Assume 'data' is your DataFrame with log fold changes and 'p_values' with corresponding p-values
data = pd.DataFrame(np.random.normal(0, 1, (1000, 5)), columns=['Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5'])
p_values = pd.DataFrame(np.random.uniform(0, 1, (1000, 5)), columns=['Cond1', 'Cond2', 'Cond3', 'Cond4', 'Cond5'])
breakpoint()

# Filtering or adjusting data based on p-values
significance_threshold = 0.05
data_filtered = data.where(p_values <= significance_threshold)

# Optional: Replace NaNs with 0 or mean or any other imputation method
data_filtered.fillna(0, inplace=True)

# Hierarchical clustering
linked = linkage(pdist(data_filtered, 'euclidean'), method='ward')
breakpoint()
plt.figure(figsize=(10, 7))
dendrogram(linked, labels=data_filtered.index.tolist(), leaf_rotation=90)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Index')
plt.ylabel('Euclidean Distance')
plt.show()
breakpoint()
