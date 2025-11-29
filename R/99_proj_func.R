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


## Function for renaming replicates in Figure2a (readability and comparability) 

strip_pattern <- function(df, pattern) {
  # pattern : pattern to remove from replicates names 
  df |>
    mutate(
      rep1 = ifelse(grepl(pattern, rep1),
                    stringr::str_remove(rep1, pattern),
                    rep1),
      rep2 = ifelse(grepl(pattern, rep2),
                    stringr::str_remove(rep2, pattern),
                    rep2)
    )
}


## Function that takes a correlation matrix and outputs the long-format dataframe

cor_mat_to_long <- function(corr_mat, rowname_col = "rep1",
                            value_col = "correlation") {
  corr_mat |>
    as.data.frame() |>
    rownames_to_column(rowname_col) |>
    pivot_longer(
      cols = -all_of(rowname_col),
      names_to = "rep2",
      values_to = value_col
    )
}
