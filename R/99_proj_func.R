#Function to calculate the percentage of NA's in a data frame
NA_check <- function(dataframe) {
  total_values = nrow(dataframe) * ncol(dataframe) 
  total_NAs = sum(is.na(dataframe))
  NA_percentage = (total_NAs/total_values) * 100 
  NA_percentage 
}


# Log2(x+1) function for tables 3 and 5 for heatmap rna and prot (2C plot)
log2p1 <- function(x) {log2(x + 1)}




# Function for calculating a mean across replicates for specific cell type  

# with this function we replaced this code:
# prot_for_plot <- clean_5 |>
#  mutate(across(adult_microglia_1:oligodendrocytes_div4_3, log2p1)) |> # applying log2
#  rowwise() |> # dplyr function that applies the next function to each row, otherwise it would be applied to columns
#  mutate( 
#    prot_microglia   = mean(c_across(starts_with("adult_microglia_")),        na.rm = TRUE), # ignoring missing values
#    prot_astro       = mean(c_across(starts_with("astrocytes_")),             na.rm = TRUE), # c_across is c(across())

mean_across_col <- function(df, prefixes) {
  # prefixes : named character vector
  # c(prot_microglia = "adult_microglia_",
  # prot_astro     = "astrocytes_")...

  # in prefix vector a value (prefix, old column name start) is stored and a name (new column name)
  for (new_name in names(prefixes)) { # here we iterate over the names
    pref <- prefixes[[new_name]] # here we retrieve the prefix string that corresponds to that name
    
    df <- df |> 
      rowwise() |> 
      mutate(mean_val_temp = mean(c_across(starts_with(pref)), na.rm = TRUE)) |> 
      ungroup()
    
    df[[new_name]] <- df$mean_val_temp
    
    # removing the temporary column
    df$mean_val_temp <- NULL
  }
  
  df
}


## Function that takes a correlation matrix and outputs the clustered long-format dataframe

cor_mat_to_long <- function(cor_mat, rowname_col = "rep1",
                            value_col = "correlation") {
  
  cor_mat |>
    as.data.frame() |>
    rownames_to_column(rowname_col) |>
    pivot_longer(
      cols = -all_of(rowname_col),
      names_to = "rep2",
      values_to = value_col
    )
}

## Function to cluster columns in a correlation matrix
cluster_cor_mat <- function(cor_mat) {
  
  # cluster using 1 - correlation
  hc <- hclust(as.dist(1 - cor_mat))
  
  # convert to long
  cor_long <- cor_mat_to_long(cor_mat)
  
  # apply cluster order
  cor_long <- cor_long |>
    mutate(
      rep1 = factor(rep1, levels = hc$labels[hc$order]),
      rep2 = factor(rep2, levels = hc$labels[hc$order])
    )
  
  return(cor_long)
}


## Function to plot rank-rank comparison with specific set of genes
plot_rank_rank_comparison <- function(
    marker_list,        # Vector of gene names to highlight
    cell_type_name,     # String for the plot title (e.g., "Astrocyte")
    use_density = FALSE,# Boolean: TRUE to use density contours, FALSE for scatter plot
    plot_limit = 5000   # Maximum rank limit for both axes (removes top artifact)
) {
  
  marker_data <- merged_ranked |>
    filter(gene %in% marker_list)
  
  p <- merged_ranked |>
    ggplot(aes(x = rna_rank, y = prot_rank)) +
    
    {if(use_density) {
      list(
        geom_point(alpha = 0.1, color = "black", size = 0.5),
        geom_density_2d(color = "darkblue")
      )
    } else {
      geom_point(alpha = 0.5, color = "black", size = 1)
    }} +
    
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
    
    geom_point(data = marker_data, color = "red", size = 3) +
    geom_text_repel(
      data = marker_data, 
      aes(label = gene),
      color = "red",
      size = 4, 
      max.overlaps = Inf, 
      segment.color = 'gray50'
    ) +
    
    scale_x_continuous(limits = c(0, plot_limit), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, plot_limit), expand = c(0, 0)) + 
    coord_equal() + 
    
    labs(
      title = paste("Rank-Rank Comparison:", cell_type_name, "RNA vs. Protein"),
      x = "RNA Rank (Low Rank = High Expression)",
      y = "Protein Rank (Low Rank = High Abundance)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0("../results/rank_", cell_type_name, ".jpg"),
    plot = p, 
    width = 18,
    height = 10, 
    dpi = 300)
  
  return(p)
}