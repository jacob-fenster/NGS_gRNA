# Load necessary libraries
library(tidyverse)
library(ggpubr)
library(gridExtra)

# Read the data from a TSV file
data <- read_tsv("/Users/jfenster/Documents/NGS_gRNA/output/20231126_exact_counts/20231126_Carbon_exact_MAGeCK.txt")

# Normalize counts by row
normalized_data <- data %>%
  mutate(across(starts_with("LB_") | starts_with("gluc_") | starts_with("ace_") | starts_with("10pCA_") | starts_with("50pCA_"), 
                ~ . / rowSums(across(starts_with("LB_") | starts_with("gluc_") | starts_with("ace_") | starts_with("10pCA_") | starts_with("50pCA_")))))

# Function to create scatter plots with linear fit and R² value
create_scatter_plot <- function(data, condition) {
  condition_R1 <- sym(paste0(condition, "_R1"))
  condition_R2 <- sym(paste0(condition, "_R2"))
  
  plot <- ggplot(data, aes(x = !!condition_R1, y = !!condition_R2)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    stat_cor(method = "pearson", label.x.npc = "right", label.y.npc = "top") +
    labs(title = paste("Scatter plot of", condition, "R1 vs R2"),
         x = paste(condition, "R1"), y = paste(condition, "R2"))
  return(plot)
}

# Plot scatter plots for all conditions
conditions <- c("LB_t0", "LB_t1", "LB_t2", "gluc_t1", "gluc_t2", "ace_t1", "ace_t2", "10pCA_t1", "10pCA_t2", "50pCA_t1", "50pCA_t2")

plots <- lapply(conditions, create_scatter_plot, data = normalized_data)

# Display plots
do.call(grid.arrange, c(plots, ncol = 2))

