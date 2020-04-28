rm(list=ls())
# for checking package availability, if not, install it
source("R/package_verification.R")

source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/read_gcmc_data.R")
source("R/read_textural_data.R")
source("R/save_train_test_data.R")
source("R/refined_bins_calc.R")
# in command prompt: Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ use_ch4
# tell it whether to save plots for poster
poster <<- TRUE # save this as a global variable: no need to pass it around
# poster plots should have big dots
if (poster){
  dot_size <<- 3
}else{
  dot_size <<- 1 
}
# define the input files
gcmc_file = "All_data/CoRE_Hexane_.1bar_298K.txt" # try the results from the CBCFC simulations
grid_file = "All_data/CH3_1A_combined_tobacco_core.rds"
previous_plot_lim <- 200
# extract the directory from the file name
data_dir <- sub("\\/.*", "", grid_file)
file_name <- sub(".*\\/", "", grid_file)
data_dir <- getwd() %>% paste0(., "/", data_dir, "/")
# find files with that prefix
# because 0.5 A is too large, we splitted it into 5 rds files
list_of_files <- list.files(data_dir, pattern = file_name)
# then check length of the list of files
if (length(list_of_files) > 1){
  for (i in 1:length(list_of_files)){
    gridi = read_rds(paste0(grid_file, as.character(i), ".rds"))
    if (i == 1) {
      grids = gridi
    } else{
      grids <- rbind(grids, gridi)
    }
  }
  hmof_h2_grid <- grids
} else{
  hmof_h2_grid <- read_rds(grid_file)
}

if (!dir.exists("Results")){
  dir.create("Results")
}
#unit_list <- c("Molec_cm3overcm3", "Mol_kg", "Mill_gram", "Cm3_gram")
unit_list <- c("Molec_cm3overcm3")
for (chosen_unit in unit_list){ 
  gcmc_data <- read_data(gcmc_file, unit_of_ads = chosen_unit, no_low_loading = FALSE)
  # what if we filter those zero loading structures out?
  gcmc_data <- gcmc_data %>% filter(Uptake > 1)
  # extract the probe name, molecule name, pressure and temperature for picture naming
  # the naming convention of rds file should be: "probe"_"size of grid"_whatever.rds
  if (str_count(grid_file, "_") > 1) {
    probe_name <- sub(".*\\/", "", grid_file) %>% sub("\\_.*", "", .)
    grid_size <- sub(".*\\/", "", grid_file) %>% 
      str_match(., paste0(probe_name, "_(.*?)_"))
    grid_size <- grid_size[2]
    grid_info <- paste(probe_name, grid_size, sep = "_")
  } else{
    grid_info <- sub(".*\\/", "", grid_file) %>% sub("\\.rds", "", .)
  }
  molecule_name <- sub(".*\\/", "", gcmc_file) %>% sub("\\_.*", "", .)
  temperature <- gcmc_data$Temp[1]
  pressure <- gcmc_data$Pres[1]
  string_to_paste <- paste(molecule_name, temperature, pressure, grid_info, sep = "_")
  # create a directory for that condition and molecule
  save_path <- paste0(string_to_paste, "/", sep = "")
  # if save for poster figures, put poster in the front
  if(poster){
    save_path <- paste0("Poster_", save_path)
  }
  save_path <- paste0("Results/", save_path)
  if (!dir.exists(save_path)){
    dir.create(save_path)
  }
  # now, we want to consider the effect of units, so add units to save_path
  unit_path <- paste0(save_path, chosen_unit, "/")
  if (!dir.exists(unit_path)){
    dir.create(unit_path)
  }
  save_path <- unit_path
  # here we can have two options, feed uniform or non-uniform bins
  #binbounds <- bounds_from_params(ch4_binspec)
  # # delete the line below, just testing Core MOFs
  # a <- read_table2(gcmc_file)
  # b <- a %>% select(-Heat_of_Ads, -Heat_fluc)
  # bb <- na.omit(b)
  # gcmc_data <- bb %>% select(ID, Temp, Pres, Molec_cm3overcm3)
  # gcmc_data <- mutate(gcmc_data, Uptake = Molec_cm3overcm3)# delete!
  gcmc_data <- mutate(gcmc_data, id=ID)
  new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
  Uniform_bins <- FALSE
  if (!Uniform_bins){
    binbounds <- automatic_bins(new_grid) # for non-uniform bins
    binbounds <- as.data.frame(binbounds)
    binbounds <- binbounds %>% filter(upper <= 0.0) # bind all repulsive bins
  } else{
    # if we want to use uniform bins
    binbounds <- data.frame()
    # 1 KJ/mol spacing
    lower <- seq(-50, -1, by = 1)
    upper <- seq(-50+1, 0, by = 1)
    binbounds <- data.frame(lower, upper)
  }
  ch4_binspec <- list("from" = head(binbounds$lower, n = 1), 
                      "to" = tail(binbounds$upper, n = 1), 
                      bounds = binbounds)
  
  
  # depending on the size of the dataset, data split is going to have two options
  if (nrow(gcmc_data) >= 2000){
    DATA_SPLIT <- 1000 # Number of data points used in training data, from setup_data.R
  } else {
    DATA_SPLIT <- ceiling(0.5 * (nrow(gcmc_data)))
  }
  
  p_ch4_sets <- partition_data_subsets(hmof_h2_grid, gcmc_data, DATA_SPLIT)
  p_ch4_vol <- gcmc_data %>%
    mutate(g.L = Uptake) %>%
    run_model_on_partitions(p_ch4_sets, ., ch4_binspec, 
                            plot_units=unit_for_plot, db_name = "tobacco")
  # calculate the plot limit: largest value in the dataset
  plot_limit <- max(gcmc_data$Uptake)
  # round to the nearest 10
  plot_limit <- round(plot_limit, digits = -1) + 10
  gg_train <- p_ch4_vol$plots$parity_training + 
    theme(axis.title.y = element_text(hjust=0.5)) + 
    xlab(paste0("GCMC capacity",  "(", unit_for_plot, ")")) + 
    ylab(paste0("Predicted capacity", "(", unit_for_plot , ")"))
  
  gg_test <- p_ch4_vol$plots$parity_testing + 
    theme(axis.title.y = element_text(hjust=0.5)) + 
    xlab(paste0("GCMC capacity",  "(", unit_for_plot, ")")) + 
    ylab(paste0("Predicted capacity", "(", unit_for_plot , ")"))
  if (poster){
    gg_train <- gg_train + theme(axis.text=element_text(size=30),
                                 axis.title=element_text(size=30,face="bold"))
    gg_test <- gg_test + theme(axis.text=element_text(size=30),
                               axis.title=element_text(size=30,face="bold")) 
  }
  save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_train.png"), sep = ""), 
            gg_train, base_width = 10, base_height = 10, dpi = 600)
  save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_test.png"), sep = ""), 
            gg_test, base_width = 10, base_height = 10, dpi = 600)
  # export data for later testing other models, one can also export data to csv files
  train_data <- export_data(p_ch4_sets$training, 
                            rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), ch4_binspec)
  test_data <- export_data(p_ch4_sets$testing, 
                           rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), ch4_binspec)
  
  rbind(train_data, test_data) %>% write_csv("latest_tob_simulation_1e6_041920.csv")
  # a histogram correlation analysis
  train_histo_bins <- train_data %>% select(-y_act)
  test_histo_bins <- test_data %>% select(-y_act)
  cor(as.matrix(t(train_histo_bins[1:2,])))
  # Run a random forest model
  rf_model <- randomForest(x = train_data %>% select(-y_act), y = train_data$y_act, ntree = 500)
  # gg <- qplot(x = model$predicted, y = model$y) + geom_abline(slope = 1, intercept = 0)
  
  # # plot error automatically if model is a randomforest
  # #plot(model)
  # # making predictions/ testing your fit, tell what variables are important
  #varplot <- varImpPlot(rf_model)
  #save_plot(paste("Results/", paste0(string_to_paste, "_Variable_Importance_Plot.png"), sep = ""),varplot, base_width = 10, base_height = 8, dpi = 600)
  # 
  make_rf_prediction_plots(condition_name = string_to_paste, 
                           plot_name = "rf", lim = plot_limit, 
                           rf_model= rf_model,  test_data = test_data)
  
  # Part below are related to topology, see Paper by Yamil and Diego
  # Read the names with topology
  topologies <- read.table("All_data/fullnames_without_tob_cleaner_py.txt") # first column is ID, second column is topology
  
  # just keep all the structural properties that are numbers
  # filter out those tobacco data with NAs
  ToBaCCo <- TRUE
  CoRE <- TRUE
  if (ToBaCCo & CoRE){
    data_1 <- na.omit(read_textual_data())
    data_2 <- na.omit(read_textual_data(option = "CoRE"))
    orig_tobacco_data <- rbind(data_1, data_2)
  }
  else if (ToBaCCo){
    orig_tobacco_data <-  na.omit(read_textual_data())
  }
  else if(CoRE){
    orig_tobacco_data <-  na.omit(read_textual_data(option = "CoRE"))
  }
  #
  structural_data <- orig_tobacco_data # %>% select(MOF.ID, vf, vsa, gsa, pld, lcd)
  
  
  train_data_with_id <- export_data(p_ch4_sets$training, 
                                    rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), 
                                    ch4_binspec, with_id = TRUE)
  test_data_with_id <- export_data(p_ch4_sets$testing, 
                                   rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), 
                                   ch4_binspec, with_id = TRUE)
  # calculate the correlation matrix for these data
  colnames(train_data_with_id)[colnames(train_data_with_id)=="id"] <- "MOF.ID"
  colnames(test_data_with_id)[colnames(test_data_with_id)=="id"] <- "MOF.ID"
  correlation_train <- correlation_of_energy_histograms(train_data_with_id)
  correlation_test <- correlation_of_energy_histograms(test_data_with_id)
  Whole_correlation <- correlation_of_energy_histograms(rbind(train_data_with_id, test_data_with_id), write_to_csv = TRUE)
  # add R-score as a predictor
  all_data_with_R <- rbind(train_data_with_id, test_data_with_id)
  correlation_all <- correlation_of_energy_histograms(all_data_with_R)
  # randomly select an ID as the reference for R-scores
  random_ID <- sample(all_data_with_R$MOF.ID, size = 1)
  all_Rs <- subset(correlation_all, rownames(correlation_all) %in% random_ID)
  df_all_Rs <- data.frame(t(all_Rs))
  colnames(df_all_Rs)[colnames(df_all_Rs)==paste0("X", random_ID)] <- "R_scores"
  all_data_with_R$Rscore <- df_all_Rs$R_scores
  training_rows <- sample(length(all_data_with_R$MOF.ID), DATA_SPLIT)
  training_data_with_R <- all_data_with_R[training_rows, ]
  testing_data_with_R <- all_data_with_R[-training_rows, ]
  # show that Rscore and uptake are related
  #Rscore_compare <- all_data_with_R %>% select(MOF.ID, y_act, Rscore)
  #Rscore_2 <- Rscore_compare[order(-Rscore_compare$Rscore),]
  
  
  rf_model_Rscore <- randomForest(x = training_data_with_R %>% select(-y_act, -MOF.ID),
                                  y = training_data_with_R$y_act, ntree = 500)
  make_rf_prediction_plots(condition_name = string_to_paste, 
                           plot_name = "rf_histogram_Rscore", lim = plot_limit, 
                           rf_model= rf_model_Rscore,  test_data = testing_data_with_R)
  
  # learning using Topology data
  topo_train_data<- merge(train_data_with_id, structural_data, by ="MOF.ID")
  topo_test_data <- merge(test_data_with_id, structural_data, by ="MOF.ID")
  
  # perform a random forest model with these structural properties
  # check if training and testing data has MOF.ID that has NAs in tobacco data
  
  rf_model_topo <- randomForest(x = topo_train_data %>% select(-y_act, -MOF.ID), 
                                y = topo_train_data$y_act, ntree = 500)
  make_rf_prediction_plots(condition_name = string_to_paste, 
                           plot_name = "rf_histogram_topology", lim = plot_limit, 
                           rf_model= rf_model_topo,  test_data = topo_test_data)
}

rbind(topo_test_data, topo_train_data) %>% write_csv("latest_tob_simulation_topo_1e6_041920.csv")


# for outliers in hexane simulations
# if (grepl("Hexane_298_10000", string_to_paste))
# {
#   # also make histograms for those large deviated points
#   # filter out those whose difference between predicted and actual is greater than 40
#   tested <- predict(rf_model_topo, topo_test_data)
#   train_pred <- predict(rf_model_topo)
#   large_deviations_train <- topo_train_data
#   large_deviations_test <- topo_test_data
#   large_deviations_train$pred <- train_pred
#   large_deviations_test$pred <- tested
#   large_deviations <- rbind(large_deviations_train, large_deviations_test)
#   small_deviations <- large_deviations %>% filter(abs(y_act - pred) < 30)
#   large_deviations <- large_deviations %>% filter(abs(y_act - pred) >= 30)
#   axis_label <- "capacity (cm\u00B3/cm\u00B3)"
#   gg_rf_test_without_outliers <- qplot(x = small_deviations$y_act, y = small_deviations$pred) +
#     geom_abline(slope = 1, intercept = 0, linetype = 2) +
#     xlab(paste0("GCMC ", axis_label)) +
#     ylab(paste0("Predicted ", axis_label))
#   # get the gcmc_data for those who has large deviations
#   large_deviation_gcmc <- gcmc_data %>% filter(ID %in% large_deviations$MOF.ID)
#   # filter the original tobacco data
#   large_deviation_tobacco <- orig_tobacco_data %>% filter(MOF.ID %in% large_deviations$MOF.ID)
#   new_string_to_paste <- paste0(string_to_paste, "_2")
#   make_topology_histograms(condition_name = new_string_to_paste, topo_data = large_deviation_tobacco)
#   # now lets distinguish between outlier: overpredict vs. underpredict
#   large_deviation_positive <- large_deviations %>% filter((y_act - pred) > 0)
#   large_deviation_negative <- large_deviations %>% filter((y_act - pred) <= 0)
# }
# large_deviation_positive %>% write_csv("outlier_under_with_pred.csv")
# large_deviation_negative %>% write_csv("outlier_over_with_pred.csv")
# # interpret the rf model
# treelist <- RF2List(rf_model_topo)
# exec <- extractRules(treelist, topo_train_data %>% select(-y_act, -MOF.ID))
# exec[1:2,]
# # get rule metric
# ruleMetric <- getRuleMetric(exec,topo_train_data %>% select(-y_act, -MOF.ID),topo_train_data$y_act)
# ruleMetric[1:2,]
# # try getTree
# tree <- getTree(rf_model_topo, 1, labelVar = TRUE)
# 
# to.dendrogram <- function(dfrep,rownum=1,height.increment=0.1){
#   
#   if(dfrep[rownum,'status'] == -1){
#     rval <- list()
#     
#     attr(rval,"members") <- 1
#     attr(rval,"height") <- 0.0
#     attr(rval,"label") <- dfrep[rownum,'prediction']
#     attr(rval,"leaf") <- TRUE
#     
#   }else{##note the change "to.dendrogram" and not "to.dendogram"
#     left <- to.dendrogram(dfrep,dfrep[rownum,'left daughter'],height.increment)
#     right <- to.dendrogram(dfrep,dfrep[rownum,'right daughter'],height.increment)
#     rval <- list(left,right)
#     
#     attr(rval,"members") <- attr(left,"members") + attr(right,"members")
#     attr(rval,"height") <- max(attr(left,"height"),attr(right,"height")) + height.increment
#     attr(rval,"leaf") <- FALSE
#     attr(rval,"edgetext") <- dfrep[rownum,'split var']
#     #To add Split Point in Dendrogram
#     #attr(rval,"edgetext") <- paste(dfrep[rownum,'split var'],"\n<",round(dfrep[rownum,'split point'], digits = 2),"=>", sep = " ")
#   }
#   
#   class(rval) <- "dendrogram"
#   
#   return(rval)
# }
# 
# d <- to.dendrogram(tree)
# str(d)
# plot(d,center=TRUE,leaflab='none',edgePar=list(t.cex=1,p.col=NA,p.lty=0))
