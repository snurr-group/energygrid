# save all the data to a csv file
# export_training_data <- function(partitioned_mod, filename) {
#   # %$% exposition pipe operator, magrittr package
#   partitioned_mod$trained_mod %$% bind_cols(y=y, orig_x) %>% write_csv(filename)
#   
# }


# Now get the testing data, using a combination of pred_grid/eval_test_grid.
# convert it to a general function, not just for test data
# can combine with other data like heat of adsorption
export_data <- function(test_grid, df_with_y_act, binspec, with_id = FALSE, align_bins="strict", writecsv = FALSE, filename = "default_data.csv") {
  grid_desc <- test_grid %>%
    stepped_hist(binspec, align_bins = align_bins) %>% 
    spread(key=bin, value=metric)
  results <- df_with_y_act %>% 
    select(id, y_act) %>%
    inner_join(grid_desc, by="id")
  if (!with_id){
    results <- results %>% select(-id)
  }
  if (writecsv){
    results %>% write_csv(filename)
  }
  results
}
# extra function for extracting a separate set of textural properties data from other training purposes
# Just toplogies
export_textural_data <- function(test_grid, df_with_y_act,text_data,filename){
  # first get those IDs
  new_df <- subset(df_with_y_act, df_with_y_act$id %in% test_grid$id)
  # join it with gcmc_data
  text_data$MOF.ID <- as.character(text_data$MOF.ID)
  colnames(text_data)[colnames(text_data) == "MOF.ID"] <- "id"
  
  results <- text_data %>% 
  select(id,VF,VSA,GSA,PLD,LCD) %>% 
  inner_join(new_df, by = "id") %>% 
  select(`Guest-Guest`,VF,VSA,GSA,PLD,LCD)
  results %>% write_csv(filename)
}