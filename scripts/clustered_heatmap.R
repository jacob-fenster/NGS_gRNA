# This is run for the gene data
library(dplyr)
library(ComplexHeatmap)
library(glue)
library(circlize)

## global parameters 
data = read.delim("/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/clustering/20240722_gene_lfc.tsv",
                  row.names=1)
## Create the data matrix where columns are growth conditions and rows are genes
data_matrix <- as.matrix(data)
condition_types <- c("LB", "glucose", "acetate", "low_pCA", "high_pCA")
# Create a data frame for the top annotation
annotation_df <- data.frame(Condition = condition_types)

# Define colors for the annotations (optional)
annotation_colors <- list(Condition = c(LB = "blue", glucose = "green", acetate = "orange", low_pCA = "purple", high_pCA = "red"))
# Define the color scale
col_fun <- colorRamp2(c(-1, -0.8, 0, 0.8, 1), c("darkblue", "#0084ff", "white", "#ffaa00", "red"))

# Define the heatmap
heatmap <- Heatmap(data_matrix,
                   name = "Log Fold Change",
                   row_title = "Genes",
                   column_title = "Conditions",
                   row_dend_reorder = TRUE,  # Reorder the dendrogram based on row clustering
                   column_dend_reorder = TRUE,  # Reorder the dendrogram based on column clustering
                   cluster_rows = TRUE,  # Cluster genes
                   cluster_columns = TRUE,  # Cluster conditions
                   top_annotation = HeatmapAnnotation(df = annotation_df,
                                                      col = annotation_colors),
                   show_row_names = FALSE,
                   show_column_names = TRUE,
                   show_row_dend = TRUE,
                   show_column_dend = TRUE,
                   col = col_fun)

# Draw the heatmap
pdf("heatmap.pdf", width = 10, height = 8)
draw(heatmap)
dev.off()

# Save the heatmap as a PNG file
png("heatmap.png", width = 1000, height = 800, res = 150)
draw(heatmap)
dev.off()
# Extract the order of rows and columns after clustering
row_order <- row_order(heatmap)
column_order <- column_order(heatmap)

# Reorder the data matrix based on the clustering results
clustered_data <- data_matrix[row_order, column_order]

# Convert the clustered data matrix to a data frame for inspection
clustered_data_df <- as.data.frame(clustered_data)

# Output the clustered data for inspection
write.table(clustered_data_df, file = "clustered_data.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

# Print the clustered data to the console for quick inspection
#print(clustered_data_df)