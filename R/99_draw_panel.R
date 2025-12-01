# Function to plot the log2 fold expression for a given cell type
# Function that draws a single panel for a given cell type
# df_sub    = data for that cell type
# df_labels = only the marker genes for labeling
# panel_title = title shown above the panel

draw_panel <- function(df_sub, df_labels, panel_title) {
  
  # Initializing the ggplot using the subset of data for the cell type (df_sub)
  # logFC is plotted on the x-axis and logLFQ on the y-axis
  ggplot(df_sub, aes(x = logFC, y = logLFQ)) +
    
    # Draw all proteins as points, colored by their assigned category (depending on cell type)
    geom_point(aes(color = category), alpha = 0.95) +
    
    # Using color palette defined earlier
    scale_color_manual(values = cols) +
    
    # Adding a vertical dashed line at x = 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    
    # Adding gene labels only for marker genes
    # df_labels contains only the selected markers for this cell type
    geom_text_repel(
      data = df_labels,
      aes(label = gene),
      size = 3
    ) +
    
    # Adding the panel title and axis labels
    labs(
      title = panel_title,
      x = "log2 fold expression (cell type vs median)",
      y = "log2 LFQ intensity"
    ) +
    
    # Apply the shared theme used across all panels
    base_theme
}
