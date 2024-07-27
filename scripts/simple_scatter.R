library(ggplot2)

# Function to create and save a scatter plot using ggplot2 with tidy evaluation
create_scatter_plot <- function(data, x_col, y_col, filename) {
  # Check if specified columns exist in the data frame
  browser()
  if (!(x_col %in% names(data)) || !(y_col %in% names(data))) {
    stop("Specified columns do not exist in the data frame")
  }
  
  # Create the scatter plot using tidy evaluation
  plot <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point() + # Adds points to the plot
    theme_minimal() + # Optional: Uses a minimal theme for the plot
    labs(x = x_col, y = y_col, title = paste("Scatter Plot of", x_col, "vs", y_col))
  
  # Save the plot to a file
  ggsave(filename, plot, width = 7, height = 5)
}

# Example usage
# Assuming 'data' is your data frame, and it has columns 'height' and 'weight'
my_data <- data.frame(height = c(160, 170, 180, 190), weight = c(60, 70, 80, 90))
create_scatter_plot(my_data, "height", "weight", "scatter_plot.png")
