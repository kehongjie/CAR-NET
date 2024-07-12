#function to filter for different expression dataset
#input: an expression matrix with column as lncRNA or gene and row as samples
#output: an expression matrix after filtering and the indx of the genes filtered out 
#arguments: type = "ncRNA"/ "gene", single_cell = TRUE/FALSE
internal_filter = function(M, fence){
  df = M
  avg_df = apply(df, 2, mean)
  ind = avg_df < fence
  return(list(df[,!ind], ind))
}





