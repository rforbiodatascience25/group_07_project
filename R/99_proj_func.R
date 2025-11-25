#Function to calculate the percentage of NA's in a data frame
NA_check <- function(dataframe) {
  total_values = nrow(dataframe) * ncol(dataframe) 
  total_NAs = sum(is.na(dataframe))
  NA_percentage = (total_NAs/total_values) * 100 
  NA_percentage 
}