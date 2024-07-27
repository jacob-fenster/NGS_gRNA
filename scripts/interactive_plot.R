library(plotly)
library(dplyr)

# Example data generation
set.seed(42) # for reproducibility
data <- data.frame(
  Gene = paste("Gene", 1:200),
  logFC = rnorm(200, 0, 2),
  pvalue = runif(200, 0, 0.05)
)

# Calculating adjusted p-values (you might use different methods like BH adjustment)
data <- data %>%
  mutate(adjPValue = p.adjust(pvalue, method = "bonferroni"),
         Significance = ifelse(adjPValue < 0.05 & abs(logFC) > 1, "Significant", "Not Significant"))

# Creating the plot
fig <- plot_ly(data, x = ~logFC, y = ~-log10(pvalue), color = ~Significance, colors = c('#FFA07A', '#20B2AA'),
               type = 'scatter', mode = 'markers', text = ~Gene, marker = list(size = 10),
               hoverinfo = 'text+x+y') %>%
  layout(title = "Volcano Plot",
         xaxis = list(title = "Log Fold Change"),
         yaxis = list(title = "-log10(p-value)"),
         hovermode = "closest")

# Adding labels on click
fig <- fig %>% event_register("plotly_click") %>%
  event_register("plotly_selected") %>%
  onRender("
    function(el, x) {
      var trace = x.data[0];
      el.on('plotly_click', function(data) {
        var pts = data.points[0];
        var newAnnotation = {
          x: pts.x,
          y: pts.y,
          xref: 'x',
          yref: 'y',
          text: pts.text,
          showarrow: true,
          arrowhead: 7,
          ax: 0,
          ay: -40
        };
        var div = document.getElementById(el.id);
        var newIndex = div.layout.annotations ? div.layout.annotations.length : 0;
        Plotly.relayout(el.id, 'annotations[' + newIndex + ']', newAnnotation);
      });
    }
  ")

# Display the plot
fig
