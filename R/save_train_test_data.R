library(magrittr)
# save all the data to a csv file
# export_training_data <- function(partitioned_mod, filename) {
#   # %$% exposition pipe operator, magrittr package
#   partitioned_mod$trained_mod %$% bind_cols(y=y, orig_x) %>% write_csv(filename)
#   
# }


# Now get the testing data, using a combination of pred_grid/eval_test_grid.
# convert it to a general function, not just for test data
# can combine with other data like heat of adsorption
export_data <- function(test_grid, df_with_y_act, filename, binspec, align_bins="strict") {
  grid_desc <- test_grid %>%
    stepped_hist_spec(binspec, align_bins = align_bins) %>% 
    spread(key=bin, value=metric)
  results <- df_with_y_act %>% 
    select(id, y_act) %>% 
    inner_join(grid_desc, by="id") %>% 
    select(-id)
  results %>% write_csv(filename)
}