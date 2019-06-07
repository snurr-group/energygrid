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
    stepped_hist(binspec, align_bins = align_bins) %>% 
    spread(key=bin, value=metric)
  results <- df_with_y_act %>% 
#    select(id, y_act) %>%
    select(id, y_act, Heat_of_Ads,`Host-Guest`,`Guest-Guest`,Qst_low_loading) %>% 
    inner_join(grid_desc, by="id") %>% 
    select(-id)
  results %>% write_csv(filename)
}
# extra function for extracting a separate set of topology data from other training purposes
# Just toplogies
export_topology_data <- function(test_grid, df_with_y_act,topo_data,filename){
  # first get those IDs
  new_df <- subset(df_with_y_act, df_with_y_act$id %in% test_grid$id)
  # join it with gcmc_data
  topo_data$MOF.ID <- as.character(topo_data$MOF.ID)
  colnames(topo_data)[colnames(topo_data) == "MOF.ID"] <- "id"
  
  results <- topo_data %>% 
  select(id,vf,vsa,gsa,pld,lcd) %>% 
  inner_join(new_df, by = "id") %>% 
  select(`Guest-Guest`,vf,vsa,gsa,pld,lcd)
  results %>% write_csv(filename)
}
