# Load necessary libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plotly))
library(htmlwidgets)

# Function to create a combined volcano plot with formatting for a ppt slide
combined_volcano_plot_ppt <- function(datasets, logFC_cols, pval_cols, labels,
                                      x_limits, y_limits, x_breaks, y_breaks) {
  # Combine datasets into a single data frame
  combined_data <- bind_rows(lapply(seq_along(datasets), function(i) {
    data <- datasets[[i]]
    logFC_col <- logFC_cols[[i]]
    pval_col <- pval_cols[[i]]
    
    # Check if the provided columns exist in the data
    if (!(logFC_col %in% colnames(data)) || !(pval_col %in% colnames(data))) {
      stop("Provided column names do not exist in the data.")
    }
    
    data %>%
      mutate(
        neg_log10_pval = -log10(!!sym(pval_col)),
        dataset = labels[[i]],
        gene = data$gene  # Assuming 'gene' is the column name for labels
      ) %>%
      select(all_of(logFC_col), neg_log10_pval, dataset, gene) %>%
      rename(logFC = !!sym(logFC_col))
  }))
  
  # Create combined volcano plot
  p <- ggplot(combined_data, aes(x = logFC, y = neg_log10_pval, text = gene)) +
    geom_point(color = "black", shape = 16, size = 1, alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Combined Volcano Plot",
      x = "Log Fold Change",
      y = "-Log10(P-value)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "black"),
      axis.text = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18, color = "black"),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1),  # Add axis lines
      legend.position = "none"  # Remove the legend
    ) +
    scale_x_continuous(limits = x_limits, breaks = seq(x_limits[1], x_limits[2], by = x_breaks)) +  # Set x-axis major tick spacing
    scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = y_breaks))  # Set y-axis major tick spacing
  
  return(p)
}

# Read the data
data <- read.delim("/Users/jfenster/Documents/NGS_gRNA/output/20231226_main_analysis/50pCA_t2-vs-50pCA_t1.gene_summary.txt")

# Define column names for log fold change and p-value for each dataset
logFC_cols <- c("pos.lfc", "neg.lfc")  # Replace with your column name
pval_cols <- c("pos.p.value", "neg.p.value")    # Replace with your column name
labels <- c("Dataset 1", "Dataset 1") # Labels for each dataset

# Filter data based on logFC
pos_data <- data %>%
  filter(!!sym(logFC_cols[1]) > 0)
neg_data <- data %>%
  filter(!!sym(logFC_cols[2]) <= 0)

# Combine datasets into a list
datasets <- list(pos_data, neg_data)

# Define axis limits
x_limits <- c(-3, 3)  # Customize as needed
y_limits <- c(0, 7)  # Customize as needed
x_breaks <- 1
y_breaks <- 1

# Create combined volcano plot
volcano_plot <- combined_volcano_plot_ppt(datasets, logFC_cols, pval_cols, 
                                          labels, x_limits, y_limits, x_breaks, y_breaks)

# Convert ggplot to plotly for interactivity
interactive_plot <- ggplotly(volcano_plot, tooltip = "text")

# Save the interactive plot as an HTML file
saveWidget(interactive_plot, "volcano_plot.html", selfcontained = TRUE)

# Display the interactive plot in the RStudio Viewer pane
viewer("volcano_plot.html")

