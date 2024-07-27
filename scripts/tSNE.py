import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

# Load your dataset
out_dir = '/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/clustering/'
data = pd.read_csv('/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/clustering/20240722_gene_lfc.tsv',
                   sep='\t', index_col=0)


data = data.drop(['lfc|LB', 'lfc|low_pCA'], axis=1)
# Standardize the data
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data)

clusters = 7
km = KMeans(n_clusters=clusters, init='k-means++', algorithm="elkan",random_state=42)  
km.fit(data_scaled)
prediction = km.predict(data_scaled)
df_km = pd.DataFrame(prediction)

# Apply t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=35)
data_tsne = tsne.fit_transform(data_scaled)
df_km['x'] = data_tsne[:, 0]
df_km['y'] = data_tsne[:, 1]
# add clusters back to dataframe
data['cluster'] = prediction
data.to_csv(f"{out_dir}cluster{clusters}.tsv", sep='\t')
# Plotting the results
color_dict = {
    0: '#FFA040',
    1: '#F94040',
    2: '#5757F9',
    3: '#B856D7',
    4: '#A0FFA0',
    5: '#40B0FF',
    6: '#F9A857'
    #7: '#40F99B'
    #8: '#FF40B5',
    #9: '#B0FF40'
}

colors = [color_dict[category] for category in prediction]
plt.figure(figsize=(10,8))
for category, color in color_dict.items():
    mask = (prediction == category)
    plt.scatter(data_tsne[:, 0][mask], data_tsne[:, 1][mask], c=color)
plt.colorbar(ticks=range(clusters))
breakpoint()
plt.show()
if False:
    plt.figure(figsize=(10, 8))
    plt.scatter(data_tsne[:, 0], data_tsne[:, 1], c='blue', marker='o')
    plt.title('t-SNE visualization')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.show()
breakpoint()